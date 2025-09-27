import pickle
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import rustworkx
from biobit.core.loc import Orientation, Strand, Interval
from biobit.core.loc.mapping import ChainMap
from biobit.toolkit import repeto
from intervaltree import IntervalTree
from joblib import Parallel, delayed, Memory
from pybedtools import BedTool

import ld
import utils
from stories.RIP import pcalling

POOL = Parallel(n_jobs=-1, verbose=100, batch_size=1)

with open(ld.INSULATORS_CACHE, 'rb') as stream:
    INSULATORS = defaultdict(lambda: defaultdict(list))
    for cfind, insulators in pickle.load(stream).items():
        for seqid, (start, end), strand in insulators:
            INSULATORS[cfind][seqid, Orientation(strand)].append((start, end))


# Parse dsRNAs and save per contig-orientation to speed up the downstream processing
def cache_dsRNA(config: ld.Config) -> tuple[
    dict[tuple[str, Orientation], Path], dict[tuple[str, Orientation], list[repeto.repeats.InvRepeat]]
]:
    with open(config.dsRNA.predicted, 'rb') as stream:
        dsRNA = pickle.load(stream)

    paths, result = {}, {}
    for (seqid, orientation), hits in dsRNA.items():
        saveto = config.dsRNA.cache / "dsRNA-filtering" / f"{seqid}_{orientation}.pkl"
        paths[seqid, orientation] = saveto

        # Load the cached data if available
        if saveto.exists():
            with open(saveto, 'rb') as stream:
                result[seqid, orientation] = pickle.load(stream)
            continue

        unique = dict()
        for _, predicted in hits:
            for rna in predicted:
                key = tuple(rna.seqranges())
                if key in unique:
                    assert unique[key] == rna, rna
                else:
                    unique[key] = rna
        unique = sorted(unique.values(), key=lambda x: x.brange().start)
        result[seqid, orientation] = unique

        saveto.parent.mkdir(parents=True, exist_ok=True)
        with open(saveto, 'wb') as stream:
            pickle.dump(unique, stream, protocol=pickle.HIGHEST_PROTOCOL)
    return paths, result


def index_peaks(cmp: pcalling.Config) -> tuple[tuple[str, str], dict[tuple[str, Orientation], IntervalTree]]:
    # Load all sample-wise peaks & make the index & group them by strand/contig
    itrees = defaultdict(lambda: IntervalTree())
    for p in BedTool(cmp.reaper.raw_peaks):
        itrees[p.chrom, p.strand].addi(p.start, p.end)
    itrees = {(contig, Orientation(strand)): v for (contig, strand), v in itrees.items()}
    return (cmp.project, cmp.ind), itrees


def optimize(
        config: ld.Config, cache: Path, seqid: str, orientation: Orientation,
        sample: pcalling.Config, smpeaks: IntervalTree, insulators: list[tuple[int, int]]
) -> tuple[
    tuple[str, str], str, Orientation, list[tuple[repeto.repeats.InvRepeat, list[int]]]
]:
    # Load all dsRNA for the current seqid-orientation
    with open(cache, 'rb') as stream:
        dsRNA = pickle.load(stream)

    # Chop dsRNA into filtering groups based on proximity
    connections = np.zeros((len(dsRNA), len(dsRNA)), dtype=bool)

    units = []
    for ind, rna in enumerate(dsRNA):
        for segment in rna.left_brange(), rna.right_brange():
            units.append((ind, segment))
    units = sorted(units, key=lambda x: x[1].start)

    for first in range(len(units)):
        find, felement = units[first]
        for second in range(first + 1, len(units)):
            sind, selement = units[second]
            if selement.start - felement.end <= config.clusters.max_distance:
                connections[find, sind] = True
            else:
                break

    # Each node is a dsRNA, each edge is a connection denoting that two elements are part of the same group
    # We need to find all connected components in this graph - each connected component is a filtering group
    graph = rustworkx.PyGraph(node_count_hint=len(dsRNA), edge_count_hint=len(connections))
    graph.add_nodes_from(range(len(dsRNA)))
    graph.add_edges_from_no_data([
        (find, sind) for find in range(len(dsRNA)) for sind in range(len(dsRNA)) if connections[find, sind]
    ])
    groups = rustworkx.connected_components(graph)

    # Load the scoring tracks
    tracks = ld.invrep_scoring.ExperimentTracks(sample).open(
        seqid, utils.assembly.seqsizes(sample.organism)[seqid], orientation.to_strand()
    )

    results = []
    for payload in groups:
        # Subsample peaks to the current group
        start, end = float('inf'), float('-inf')
        for ind in payload:
            brange = dsRNA[ind].brange()
            start, end = min(start, brange.start), max(end, brange.end)

        # Create an index of all enrichment regions overlapping this bounding range
        pindex = IntervalTree()
        for h in smpeaks.overlap(start, end):
            pindex.addi(h.begin, h.end)

        # Chop dsRNA into pieces overlapping with peaks
        pieces = {}
        for ind in payload:
            rna = dsRNA[ind]
            lbrange, rbrange = rna.left_brange(), rna.right_brange()
            lov, rov = pindex.overlap(lbrange.start, lbrange.end), pindex.overlap(rbrange.start, rbrange.end)
            if len(lov) == 0 or len(rov) == 0:
                continue

            # Map each peak to dsRNA coordinates
            lchain, rchain, length = [], [], 0
            for segment in rna.segments:
                lchain.append(segment.left)
                rchain.append(segment.right)
                length += segment.left.len()
            lchain, rchain = ChainMap(Interval.merge(lchain)), ChainMap(Interval.merge(rchain))

            lmapped = []
            for x in lov:
                if (mapped := lchain.map_interval((x.begin, x.end))) is not None:
                    lmapped.append(mapped)

            rmapped = []
            for x in rov:
                if (mapped := rchain.map_interval((x.begin, x.end))) is not None:
                    rmapped.append(Interval(length - mapped.end, length - mapped.start))

            # Intersect to get only segments where both arms overlap with peaks
            # overlaps = Interval.overlap(lmapped, rmapped)

            # Slightly extend and merge to ignore small gaps and cover larger dsRNA segments
            overlaps = [
                Interval(max(0, x.start - 128), min(length, x.end + 128)) for x in lmapped + rmapped
            ]
            overlaps = Interval.merge(overlaps)
            overlaps = sorted(overlaps)

            # Select and crop dsRNA segments by each overlap
            for mapped in ld.invrep_scoring.from_dsRNA_coordinates_to_global(rna, overlaps):
                key = tuple(mapped.seqranges())
                if key in pieces:
                    assert pieces[key]['dsRNA'] == mapped
                    pieces[key]['origin'].append(ind)
                else:
                    pieces[key] = {'dsRNA': mapped, 'origin': [ind]}

        if not pieces:
            continue

        # Fetch scores for each dsRNA segment
        start = min(x[0].start for x in pieces.keys())
        end = max(x[-1].end for x in pieces.keys())
        scoring = tracks.score(start, end, insulators, pindex)

        # We can include in the solution all segments with a score >= 10% of the max observed score
        # (simple heuristic to speed up the process)
        unresolved = [(scoring.score(ir['dsRNA']), ir['dsRNA']) for ir in pieces.values()]
        maxscore = max(x[0] for x in unresolved)
        unresolved = [(score, ir) for score, ir in unresolved if score >= 0.1 * maxscore]

        # # Greedy algorithm based on scoring individual segments
        cache, solution = [], []
        while unresolved:
            # Score each segment
            cache.clear()
            for _, ir in unresolved:
                score = scoring.score(ir)

                if score >= 1e-32 and score >= 0.1 * maxscore:
                    cache.append((score, ir))

            unresolved, cache = cache, unresolved
            unresolved.sort(key=lambda x: x[0])

            # We can include in the solution all top-scored & non-overlapping segments
            # Overlapping segments must be re-scored in the next iteration
            resolved, left = [], []
            cache.clear()
            while unresolved:
                score, candidate = unresolved.pop()
                lbrange, rbrange = candidate.left_brange(), candidate.right_brange()

                conflict = Interval.overlap([lbrange, rbrange], cache)
                if len(conflict) > 0:
                    left.append([score, candidate])
                else:
                    resolved.append([score, candidate])
                    cache.extend([lbrange, rbrange])

            # Resolve the selected segments and "consume" their signal
            for _, ir in resolved:
                scoring.resolve(ir)

            solution.extend(resolved)
            unresolved = left

        # Trim the solution
        solution.sort(reverse=True, key=lambda x: x[0])

        total = sum(x[0] for x in solution)
        desired = 0.9 * total
        while solution and total > desired:
            score, _ = solution[-1]
            if total - score >= desired:
                total -= score
                solution.pop()
            else:
                break

        # Solution is a list of (ir, dsRNA indices) pairs
        results.extend([
            (ir, pieces[tuple(ir.seqranges())]['origin']) for _, ir in solution]
        )
        del scoring, score, solution, unresolved, cache, pieces

    return (sample.project, sample.ind), seqid, orientation, results


for config in ld.Config.load():
    if not config.dsRNA.predicted.exists():
        print(f"Skipping {config.ind}")
        continue

    dscache, dsRNA = cache_dsRNA(config)
    assert len(dsRNA) == len(dscache)

    # Launch the processing
    memory = Memory(location=config.dsRNA.cache, verbose=0)

    # Parse sample-wise peaks
    _job = memory.cache(index_peaks, verbose=False)
    peaks = POOL(delayed(_job)(cmp) for cmp in config.comparisons)
    peaks = dict(peaks)

    print(f"Total targets: {len(config.comparisons) * len(dsRNA)}")
    # _job = memory.cache(optimize, verbose=False, ignore=['smpeaks'])
    _job = optimize

    optimized = POOL(
        delayed(_job)(
            config, dscache[seqid, orientation], seqid, orientation,
            smpl, peaks[smpl.project, smpl.ind][seqid, orientation],
            INSULATORS[config.ind][seqid, orientation]
        )
        for smpl in config.comparisons
        for seqid, orientation in dscache.keys()
        if (seqid, orientation) in peaks[smpl.project, smpl.ind]
    )

    # Aggregate the information about dsRNA segments support across all samples
    records = {}
    for _, seqid, orientation, solution in optimized:
        for ir, inds in solution:
            for ind in inds:
                rna = dsRNA[seqid, orientation][ind]
                key = (seqid, orientation, *rna.seqranges())
                if key not in records:
                    records[key] = {'dsRNA': rna, 'solutions': []}
                records[key]['solutions'].append(ir)

    # Select segments that are supported by at least X samples
    before, after, result = sum(len(x) for x in dsRNA.values()), 0, defaultdict(list)
    for (seqid, orientation, *_), data in records.items():
        solution, tag = ld.invrep_scoring.filter_segments(
            data['dsRNA'], data['solutions'], config.dsRNA.min_replication
        )
        if solution:
            for invrep, t in zip(solution, tag):
                result[seqid, orientation].append((invrep, t))
            after += 1

    # Save the results
    track, pkl = [], defaultdict(list)
    for (contig, strand), irs in result.items():
        for ir, tag in irs:
            pkl[contig, Strand(str(strand))].append(ir)
            track.append(ir.to_bed12(contig=contig, strand=str(strand), name=str(tag)) + "\n")

    print(f"{config.ind}: dsRNA passed filtering: {after} / {before} ({after / before:%})")

    config.dsRNA.filtered.parent.mkdir(parents=True, exist_ok=True)
    with open(config.dsRNA.filtered, 'wb') as stream:
        pickle.dump(pkl, stream, protocol=pickle.HIGHEST_PROTOCOL)

    with tempfile.NamedTemporaryFile() as tmp:
        with open(tmp.name, 'w') as stream:
            stream.writelines(track)
        utils.bed.tbindex(BedTool(tmp.name).sort(), config.dsRNA.filtered.with_suffix(".bed.gz"))
