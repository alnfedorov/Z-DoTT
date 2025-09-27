import pickle
from collections import defaultdict
from typing import Any, Literal

import intervaltree
import pandas as pd
import pybedtools
from biobit.core.loc import Interval, Orientation
from biobit.toolkit import countit
from biobit.toolkit.repeto.repeats import InvRepeat
from intervaltree import IntervalTree
from joblib import Parallel, delayed
from pybedtools import BedTool

import utils
from stories import annotation, A2I
from stories.RIP.clustering.ld import config as clustering

BED_DTYPES = {0: str, 1: int, 2: int, 3: str, 4: str, 5: str}


def build_universe(config: clustering.Config) -> pd.DataFrame:
    def load_peaks(cmp: clustering.pcalling.Config):
        segments = defaultdict(list)
        for p in pybedtools.BedTool(cmp.reaper.filtered_peaks):
            segments[p.chrom, p.strand].append(Interval(p.start, p.end))
        segments = {(contig, Orientation(strand)): sgs for (contig, strand), sgs in segments.items()}
        return cmp.ind, segments

    # Load all peaks
    _peaks = Parallel(n_jobs=-1)(delayed(load_peaks)(cmp) for cmp in config.comparisons)

    # Add curated regions
    if config.peaks.curated_include.exists():
        curated = defaultdict(list)
        for p in BedTool(config.peaks.curated_include):
            curated[p.chrom, Orientation.Forward].append(Interval(p.start, p.end))
            curated[p.chrom, Orientation.Reverse].append(Interval(p.start, p.end))
        _peaks.append(('Curated', dict(curated)))
    else:
        print(f"WARNING: No curated regions found for the {config.ind} dataset")

    # Resolve into individual "pieces"
    peaks = dict(_peaks)
    resolved = []

    def resolution(contig: str, orientation: Orientation, start: int, end: int, comparisons: set[Any]) -> list[Any]:
        if comparisons:
            resolved.append({"contig": contig, "orientation": orientation, "start": start, "end": end})
        return []

    countit.utils.resolve_annotation(peaks, resolution)
    return pd.DataFrame(resolved)


def replication(data: pd.DataFrame, config: clustering.Config, which: Literal['raw', 'filtered']):
    def overlap_with_peaks(cmp: clustering.pcalling.Config) -> dict[str, int]:
        index = defaultdict(intervaltree.IntervalTree)
        path = {'raw': cmp.reaper.raw_peaks, 'filtered': cmp.reaper.filtered_peaks}[which]
        bed = pd.read_csv(path, sep='\t', header=None, dtype=BED_DTYPES)
        for contig, start, end, _, _, strand in bed.itertuples(index=False, name=None):
            index[contig, strand].addi(start, end)
        index = {(contig, Orientation(strand)): itree for (contig, strand), itree in index.items()}

        counts = defaultdict(int)
        for contig, orientation, start, end in data[
            ['contig', 'orientation', 'start', 'end']
        ].itertuples(index=False, name=None):
            if (contig, orientation) in index and index[contig, orientation].overlaps(start, end):
                counts[contig, orientation, start, end] += 1

        return counts

    counts = Parallel(n_jobs=-1)(delayed(overlap_with_peaks)(cmp) for cmp in config.comparisons)
    data['Replication'] = [
        sum(count.get(key, 0) for count in counts)
        for key in data[['contig', 'orientation', 'start', 'end']].itertuples(index=False, name=None)
    ]
    return data


def derive_covered_regions(config: clustering.Config):
    def _derive_covered_regions(path):
        bed = pd.read_csv(path, sep='\t', header=None, dtype=BED_DTYPES)
        parsed = defaultdict(list)
        for contig, start, end, _, _, strand in bed.itertuples(index=False, name=None):
            parsed[contig, strand].append(Interval(start, end))
        return dict(parsed)

    # Calculate covered regions
    _covered = Parallel(n_jobs=-1)(
        delayed(_derive_covered_regions)(path)
        for cmp in config.comparisons
        for path in [cmp.reaper.control, cmp.reaper.signal]
    )
    allkeys = set(k for d in _covered for k in d)

    covered = []
    covgroups = defaultdict(list)
    for seqid, strand in allkeys:
        allsegments = [segment for parsed in _covered for segment in parsed.get((seqid, strand), [])]
        for segment in Interval.merge_within(allsegments, distance=config.peaks.coverage_gaps_tolerance):
            covered.append(pybedtools.Interval(seqid, segment.start, segment.end, strand=str(strand)))
            covgroups[seqid, strand].append(segment)

    # Calculate all not covered regions
    not_covered = []
    seqsizes = utils.assembly.seqsizes(config.organism)
    for (seqid, strand), cov in covgroups.items():
        cov = sorted(cov, key=lambda x: x.start)
        if cov[0].start > 0:
            not_covered.append(pybedtools.Interval(seqid, 0, cov[0].start, strand=str(strand)))
        for prv, nxt in zip(cov[:-1], cov[1:]):
            not_covered.append(pybedtools.Interval(seqid, prv.end, nxt.start, strand=str(strand)))
        if cov[-1].end < seqsizes[seqid]:
            not_covered.append(pybedtools.Interval(seqid, cov[-1].end, seqsizes[seqid], strand=str(strand)))

    # Save the results
    config.peaks.covered.parent.mkdir(parents=True, exist_ok=True)
    BedTool(covered).sort().saveas(config.peaks.covered)

    config.peaks.not_covered.parent.mkdir(parents=True, exist_ok=True)
    BedTool(not_covered).sort().saveas(config.peaks.not_covered)


def editing_sites(data: pd.DataFrame, config: clustering.Config):
    host = utils.assembly.get(organism=config.host)
    allsites = defaultdict(set)

    # Add putative (rediportal) and all detected editing sites
    paths = [host.rediportal, A2I.tracks.all_passed / f"{host.name}.bed.gz"]
    for path in paths:
        assert path.exists(), path
        sites = pd.read_csv(path, sep='\t', header=None, dtype=BED_DTYPES)
        for contig, start, end, _, _, strand in sites.itertuples(index=False, name=None):
            allsites[contig, Orientation(strand)].add((start, end))

    # Create an index of editing sites
    itrees = defaultdict(lambda: intervaltree.IntervalTree())
    for (contig, orient), sites in allsites.items():
        for start, end in sites:
            itrees[contig, orient].addi(start, end)

    # Annotate peaks
    annotation = []
    for contig, orientation, start, end in data[
        ['contig', 'orientation', 'start', 'end']
    ].itertuples(index=False, name=None):
        total = 0
        for _ in itrees[contig, orientation].overlap(start, end):
            total += 1
        annotation.append(total)
    return {"Editing sites": annotation}


def repeats(data: pd.DataFrame, config: clustering.Config):
    host = utils.assembly.get(organism=config.host)
    host_contigs = set(host.seqid.sizes())

    # Construct the Repeat Masker Index
    itrees = defaultdict(lambda: intervaltree.IntervalTree())
    repmasker = pd.read_csv(host.repmasker, sep='\t', header=None, dtype=BED_DTYPES)
    for contig, start, end, name, _, _ in repmasker.itertuples(index=False, name=None):
        itrees[contig, '.'].addi(start, end, data=name)
    itrees = dict(itrees)

    # Annotate peaks
    annotation = []
    for contig, start, end in data[['contig', 'start', 'end']].itertuples(index=False, name=None):
        if contig not in host_contigs:
            annotation.append("Not host")
            continue

        weights = defaultdict(int)
        for it in itrees[contig, '.'].overlap(start, end):  # type: intervaltree.Interval
            length = min(it.end, end) - max(it.begin, start)
            assert length > 0
            weights[it.data] += length
        weights["Repeat-free"] = (end - start) - sum(weights.values())

        annotation.append(max(weights.items(), key=lambda x: x[1])[0])
    return {"RepeatMasker": annotation}


def genomic_regions(data: pd.DataFrame, config: clustering.Config):
    seqsizes = utils.assembly.seqsizes([config.host])

    itrees = defaultdict(IntervalTree)
    for name, locs in annotation.load.resolved_annotation(config.assembly).items():
        for (seqid, orient), intervals in locs.items():
            for interval in intervals:
                itrees[seqid, orient].addi(interval.start, interval.end, data=name)
    itrees = dict(itrees)

    gregions = []
    noannotation = set()
    for contig, orient, start, end in data[
        ['contig', 'orientation', 'start', 'end']
    ].itertuples(index=False, name=None):
        if contig not in seqsizes:
            gregions.append("Not host")
            continue
        elif (key := (contig, orient)) not in itrees:
            if key not in noannotation:
                noannotation.add((contig, orient))
                print(f"WARNING: No annotations for {contig}:{orient}")
            gregions.append(None)
            continue

        hits = {it.data for it in itrees[contig, orient].overlap(start, end)}
        has_intergenic, has_intronic = "intergenic" in hits, "intron" in hits
        hits -= {"intergenic", "intron"}

        if "CDS" in hits:
            gregions.append("cds")
        elif hits:
            gregions.append("exon")
        elif has_intronic:
            gregions.append("intron")
        else:
            assert has_intergenic
            gregions.append("intergenic")

    return {"Region": gregions}


def curated_regions(data: pd.DataFrame, config: clustering.Config):
    itrees = defaultdict(lambda: intervaltree.IntervalTree())
    for bed, key in (config.peaks.curated_exclude, "excluded"), (config.peaks.curated_include, "included"):
        if not bed.exists():
            continue

        for it in BedTool(bed):
            itrees[it.chrom, '.'].addi(it.start, it.end, data=key)

    excluded, included = [], []
    for contig, start, end in data[['contig', 'start', 'end']].itertuples(index=False, name=None):
        hits = {it.data for it in itrees[contig, '.'].overlap(start, end)}

        excluded.append('excluded' in hits)
        included.append('included' in hits)

    return {"Curated[excluded]": excluded, "Curated[included]": included}


def closest_neighbor_and_merged_length(data: pd.DataFrame):
    def job(contig, orientation, subdf):
        subdf = subdf.sort_values(by='start')

        # Merge touching peaks
        merged = []
        coordinates = []
        for start, end in subdf[['start', 'end']].itertuples(index=False, name=None):
            if merged and merged[-1][1] == start:
                merged[-1] = (merged[-1][0], end)
                coordinates[-1].append((start, end))
            else:
                merged.append((start, end))
                coordinates.append([(start, end)])

        # Save lengths
        length = {}
        for coord, (start, end) in zip(coordinates, merged):
            for c in coord:
                length[c] = end - start

        # Find distances
        left, right = {}, {}
        distances = [nxt - prev for (nxt, _), (_, prev) in zip(merged[1:], merged[:-1])]
        distances = [-1] + distances + [-1]
        for ind, coords in enumerate(coordinates):
            for c in coords:
                assert c not in left and c not in right
                left[c] = distances[ind]
                right[c] = distances[ind + 1]
        return (contig, orientation), {"Merged length": length, "Neighbor[left]": left, "Neighbor[right]": right}

    _distances = Parallel(n_jobs=-1)(
        delayed(job)(contig, orientation, subdf) for (contig, orientation), subdf in
        data.groupby(['contig', 'orientation'])
    )
    distances = dict(_distances)

    results = {"Merged length": [], "Neighbor[left]": [], "Neighbor[right]": []}
    for contig, orientation, start, end in data[
        ['contig', 'orientation', 'start', 'end']
    ].itertuples(index=False, name=None):
        indexed = distances[contig, orientation]
        for k, v in indexed.items():
            results[k].append(v[(start, end)])

    for k, v in results.items():
        data[k] = v
    return data


def dsRNA(data: pd.DataFrame, config: clustering.Config):
    # Load filtered dsRNAs
    with open(config.dsRNA.filtered, 'rb') as stream:
        dsRNAs = pickle.load(stream)

    def job(contig: str, orientation: Orientation, repeats: list[InvRepeat], peaks: list[Interval]):
        # Index repeats
        tree = intervaltree.IntervalTree()
        for r in repeats:
            for block in r.seqranges():
                tree.addi(block.start, block.end)

        records = {}
        # Calculate the overlap for individual peaks
        for peak in peaks:
            overlaps: set[intervaltree.Interval] = tree.overlap(peak.start, peak.end)
            if overlaps:
                records[peak] = max(x.overlap_size(peak.start, peak.end) for x in overlaps)
        return (contig, orientation), records

    # Calculate overlaps
    data['segment'] = [Interval(start, end) for start, end in data[['start', 'end']].itertuples(index=False, name=None)]
    _results = Parallel(n_jobs=-1)(
        delayed(job)(contig, orientation, dsRNAs[contig, orientation.to_strand()], subdf['segment'])
        for (contig, orientation), subdf in data.groupby(['contig', 'orientation'])
        if (contig, orientation.to_strand()) in dsRNAs
    )
    results = dict(_results)

    data['dsRNA[Individual]'] = [
        results.get((contig, orientation), {}).get(segment, 0)
        for contig, orientation, segment in
        data[['contig', 'orientation', 'segment']].itertuples(index=False, name=None)
    ]
    data = data.drop(columns=['segment'])
    return data
