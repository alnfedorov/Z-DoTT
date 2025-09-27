import pickle
from bisect import bisect_left
from collections import defaultdict
from itertools import chain

import numpy as np
import rustworkx
from biobit.core.loc import Interval, Orientation, Strand
from biobit.toolkit.repeto.repeats import InvRepeat
from intervaltree import IntervalTree
from joblib import Parallel, delayed
from pybedtools import BedTool, Interval as BedInterval

import ld
import utils
from utils.repeto import Partition

POOL = Parallel(n_jobs=-1, verbose=100)


def group_all(
        contig: str, strand: Strand, dsRNA: list[InvRepeat], peaks: list[Interval],
        transcripts: IntervalTree | None, curated: IntervalTree | None, insulators: list[int] | None, max_distance: int
) -> tuple[str, Strand, list[list[InvRepeat | Interval]]]:
    # Keep only dsRNAs that overlap with peaks
    index = IntervalTree.from_tuples([(x.start, x.end) for x in peaks])
    dsRNA = [x for x in dsRNA
             if index.overlaps(x.brange().start, x.brange().end) and
             any(index.overlaps(block.start, block.end) for block in x.seqranges())]

    # Try to map peaks and dsRNAs to manually curated regions first
    # Edge case:
    # * dsRNA is partially inside the curated region -> ignore them
    # * peak is partially inside the curated region -> raise an error
    groups = []
    if curated is not None:
        _dsRNA, _peaks = [], []
        overlaps = defaultdict(list)
        for peak in peaks:
            overlap = curated.overlap(peak.start, peak.end)
            if not overlap:
                _peaks.append(peak)
                continue
            assert len(overlap) == 1, (contig, strand, peak, overlap)
            ov = overlap.pop()
            assert ov.overlap_size(peak.start, peak.end) == peak.len(), (contig, strand, peak)
            overlaps[ov.data].append(peak)
        for rna in dsRNA:
            brange = rna.brange()
            lb, rb = rna.left_brange(), rna.right_brange()
            overlap = []
            for ov in curated.overlap(brange.start, brange.end):
                if ov.overlap_size(lb.start, lb.end) == lb.len() and ov.overlap_size(rb.start, rb.end) == rb.len():
                    overlap.append(ov)
            if not overlap:
                _dsRNA.append(rna)
                continue
            else:
                assert len(overlap) == 1, (contig, strand, brange)
                ov = overlap.pop()
                overlaps[ov.data].append(rna)

        dsRNA, peaks = _dsRNA, _peaks
        groups.extend(overlaps.values())

    # Enumerate all elements
    index = []
    all_elements: list[tuple[int, Interval]] = []
    for p in peaks:
        all_elements.append((len(index), p))
        index.append(p)

    for rna in dsRNA:
        for block in rna.seqranges():
            all_elements.append((len(index), block))
        index.append(rna)

    # Derive all possible connections
    connections = np.zeros((len(index), len(index)), dtype=bool)

    # * Elements that are close in genomic coordinates
    all_elements = sorted(all_elements, key=lambda x: x[1].start)
    insulators = sorted(set(insulators)) if insulators is not None else []
    closest_insulator = [bisect_left(insulators, x[1].start) for x in all_elements]
    for first in range(len(all_elements)):
        find, felement = all_elements[first]
        finsulator = bisect_left(insulators, felement.end)
        for second in range(first + 1, len(all_elements)):
            sind, selement = all_elements[second]
            sinsulator = closest_insulator[second]
            if selement.start - felement.end <= max_distance:
                # Only connect if they are in the same insulated region
                connections[find, sind] = sinsulator == finsulator
            else:
                break

    # * Elements that are close in transcriptomic coordinates
    if transcripts is not None:
        rna_overlap = defaultdict(list)
        for ind, element in all_elements:
            for rna in transcripts.overlap(element.start, element.end):
                rna_overlap[rna.data].append((ind, element))

        for rna, overlaps in rna_overlap.items():
            mapped = []
            for ind, element in overlaps:
                if rna.map(element) is not None:
                    mapped.append((ind, element))

            mapped = list({x[0] for x in mapped})
            for first in range(len(mapped)):
                find = mapped[first]
                for second in range(first + 1, len(mapped)):
                    sind = mapped[second]
                    connections[find, sind] = True

        # mapped = sorted(mapped, key=lambda x: x[1].start)
        # for first in range(len(mapped)):
        #     find, felement = mapped[first]
        #     for second in range(first + 1, len(mapped)):
        #         sind, selement = mapped[second]
        #         if selement.start - felement.end <= max_distance:
        #             connections[find, sind] = True
        #         else:
        #             break

    # Each node is a single element, each edge is a connection denoting that two elements are part of the same partition
    # We need to find all connected components in this graph - each connected component is a partition
    print(f"[{contig}:{strand}]Connections: {len(connections)}")
    graph = rustworkx.PyGraph(node_count_hint=len(index), edge_count_hint=len(connections))
    graph.add_nodes_from(range(len(index)))
    graph.add_edges_from_no_data([
        (find, sind) for find in range(len(index)) for sind in range(len(index)) if connections[find, sind]
    ])

    components = rustworkx.connected_components(graph)
    groups.extend([[index[ind] for ind in component] for component in components])
    return contig, strand, groups


def posprocess_groups(
        contig: str, strand: Strand, group: list[InvRepeat | Interval]
) -> tuple[str, Orientation, list[InvRepeat], list[Interval]]:
    # Split dsRNA/peaks
    dsRNA, peaks = [], []
    for item in group:
        if isinstance(item, InvRepeat):
            dsRNA.append(item)
        else:
            assert isinstance(item, Interval)
            peaks.append(item)

    # Keep only peaks that are not explained by predicted dsRNAs
    # = drop peaks that are fully covered by dsRNAs and its inner loops/bulges
    if peaks and dsRNA:
        peaks = sorted(peaks, key=lambda x: x.start)
        blocks = chain(*[[x.left_brange(), x.right_brange()] for x in dsRNA])
        blocks = sorted(blocks, key=lambda x: x.start)
        nonredundant = BedTool([BedInterval(contig, p.start, p.end) for p in peaks]) \
            .subtract(BedTool([BedInterval(contig, rna.start, rna.end) for rna in blocks]))

        peaks = [Interval(x.start, x.end) for x in nonredundant]

    return contig, strand.to_orientation(), dsRNA, peaks


with open(ld.INSULATORS_CACHE, 'rb') as stream:
    INSULATORS = pickle.load(stream)

for config in ld.Config.load():
    # Load peaks and group by contig/strand
    peaks = defaultdict(list)
    for p in BedTool(config.peaks.filtered):
        peaks[p.chrom, p.strand].append(Interval(p.start, p.end))
    peaks = {(contig, Strand(strand)): p for (contig, strand), p in peaks.items()}

    # Load dsRNAs and select only those overlapping with peaks
    with open(config.dsRNA.filtered, 'rb') as stream:
        dsRNA = pickle.load(stream)

    # Load curated groups
    curated = defaultdict(lambda: IntervalTree())
    if config.clusters.curated.exists():
        for ind, p in enumerate(BedTool(config.clusters.curated)):
            curated[p.chrom].addi(p.start, p.end, ind)

    # Prepare insulators
    insulators = defaultdict(list)
    for contig, (start, end), strand in INSULATORS[config.ind]:
        insulators[contig, Strand(strand)].append(start)
        insulators[contig, Strand(strand)].append(end)

    # Load well-annotated transcripts
    host = utils.assembly.get(organism=config.host)
    transcripts = ld.transcripta.parse(host)

    pre_groups = POOL(
        delayed(group_all)(
            seqid, stnd, dsRNA.get((seqid, stnd), []), peaks.get((seqid, stnd), []),
            transcripts.get((seqid, stnd)), curated.get(seqid), insulators.get((seqid, stnd)),
            config.clusters.max_distance
        )
        for seqid, stnd in peaks.keys()
    )
    pre_groups = POOL(
        delayed(posprocess_groups)(contig, strand, g)
        for contig, strand, groups in pre_groups
        for g in groups
    )

    # Enumerate & finalize groups
    groups = []
    for contig, orient, invrep, peaks in pre_groups:
        groups.append(Partition(f"P{len(groups)}", contig, orient, invrep, peaks))

    config.clusters.results.parent.mkdir(parents=True, exist_ok=True)
    with open(config.clusters.results, 'wb') as stream:
        pickle.dump(groups, stream, protocol=pickle.HIGHEST_PROTOCOL)

    track = []
    for p in groups:
        blocks = list(chain(*[item.seqranges() for item in p.invrep], p.peaks))
        blocks = Interval.merge(blocks)
        track.append(utils.bed.blocks.make(p.contig, blocks, p.orientation, p.ind))
    utils.bed.tbindex(BedTool(track).sort(), config.clusters.results.with_suffix(".bed.gz"))
