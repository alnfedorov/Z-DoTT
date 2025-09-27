import pickle
import tempfile
from collections import defaultdict
from itertools import batched
from multiprocessing import cpu_count
from typing import Any

import pybedtools
from biobit.core.loc import Orientation
from biobit.toolkit import repeto
from joblib import Parallel, delayed, Memory
from pybedtools import BedTool

import ld
import utils.repeto


# Run repeto to predict putative dsRNAs
def job(_: Any, batch: tuple[utils.repeto.RepetoGroup, ...], config: ld.Config):
    results = []
    contigs = utils.assembly.get(organism=config.host).seqid.sizes()

    for group in batch:
        if sum(x.len() for x in group.rois) < config.dsRNA.min_roi_overlap:
            # print(f"Group {group.contig}:{group.segment} has insufficient length")
            results.append((group, [], []))
            continue

        ctglen = contigs[group.contig]
        start, end = (max(0, group.bsegment.start - config.dsRNA.offset),
                      min(ctglen, group.bsegment.end + config.dsRNA.offset))
        # print(f"{group.contig}:{start}-{end} [{group.orientation}]")

        # Fetch the sequence from the genome
        seq: str = utils.fasta.sequence(
            utils.assembly.get(organism=config.host).fasta, group.contig, start, end, strand=str(group.orientation)
        ).upper().replace("T", "U")

        if group.orientation == Orientation.Reverse:
            seq = seq[::-1]

        # Predict all possible fold-back dsRNAs overlapping the peaks
        rois = [(r.start - start, r.end - start) for r in group.rois]
        filt = repeto.predict.Filter() \
            .set_min_score(config.dsRNA.min_score) \
            .set_min_matches(config.dsRNA.min_length, 1) \
            .set_min_roi_overlap(config.dsRNA.min_roi_overlap, 1) \
            .set_rois(rois)

        irs, scores = repeto.predict.run(seq.encode("ASCII"), filt, repeto.predict.Scoring())

        # Map coordinates back to the genome
        for ir in irs:
            ir.shift(start)

        # print(f"{group.contig}:{start}-{end} [{group.orientation}] -> {len(irs)}")
        results.append((group, irs, scores))

    return results


with open(ld.INSULATORS_CACHE, 'rb') as stream:
    INSULATORS = pickle.load(stream)

for config in ld.Config.load():
    print(f"Processing {config.ind}")
    config.dsRNA.cache.mkdir(parents=True, exist_ok=True)
    memory = Memory(location=config.dsRNA.cache, verbose=0)
    cached = memory.cache(job, ignore=['batch', 'config'], verbose=False)

    # Group all prefiltered peaks for dsRNA prediction
    allpeaks = [(p.chrom, (p.start, p.end), p.strand) for p in BedTool(config.peaks.prefiltered)]
    insulators = INSULATORS[config.ind]
    connectors = []

    # Load curated "connectors" for dsRNAs rev-comp lookup
    path = config.root / "curated" / "connectors.bed"
    if path.exists():
        for p in BedTool(path):
            connectors.append((p.chrom, (p.start, p.end), '+'))
            connectors.append((p.chrom, (p.start, p.end), '-'))

    # Make repeto groups
    groups = utils.repeto.group(allpeaks, insulators, connectors, maxdist=config.dsRNA.max_distance)

    # Save derived repeto groups (for visualization & debug)
    bed = [pybedtools.Interval(x.contig, x.bsegment.start, x.bsegment.end, strand=str(x.orientation)) for x in groups]
    BedTool(bed).sort().saveas(config.dsRNA.repeto)

    # Sort is required to guarantee reproducibility of the hashing
    groups = sorted(groups, key=lambda x: (x.bsegment.len(), x.contig, x.orientation, x.bsegment, x.rois))
    print(f'Total targets: {len(groups)}')

    # Print exceptionally long rois for manual inspection
    print("Longest targets:")
    for p in groups[-50:]:
        print(f"\t{p.contig}:{p.bsegment}\t{p.bsegment.len()}\tN={len(p.rois)}")
    total = len(groups)

    # Assuming that there are 1.5TB RAM and 96 cores available
    allresults = []
    for minlen, maxlen, cores, n in [
        (0, 50_000, 96, 64),
        (50_000, 75_000, 48, 1),
        (75_000, 100_000, 32, 1),
        (100_000, float("inf"), 12, 1),
    ]:
        selected, left = [], []
        for group in groups:
            if minlen <= group.bsegment.len() < maxlen:
                key = (config.ind, group.contig, group.orientation,
                       group.bsegment.start, group.bsegment.end, len(group.rois))
                selected.append((key, group))
            else:
                left.append(group)
        groups = left

        print(f"Batch length: {len(selected)}")
        selected = sorted(selected, key=lambda x: x[0])

        batch = Parallel(verbose=100, n_jobs=min(cpu_count(), cores), batch_size=1)(
            delayed(cached)(tuple(sorted(x[0] for x in batch)), tuple(x[1] for x in batch), config)
            for batch in batched(selected, n=n)
        )
        allresults.extend((group, irs, scores) for b in batch for group, irs, scores in b)

    if total != len(allresults):
        print(f"Expected {total} groups, got {len(allresults)}")
    # assert len(allresults) == len(groups)

    # Save all predicted dsRNA as BED/pkl
    track, pkl = set(), defaultdict(list)
    for group, irs, scores in allresults:
        pkl[(group.contig, group.orientation)].append((group, irs))
        for ir, score in zip(irs, scores):
            track.add(ir.to_bed12(contig=group.contig, strand=str(group.orientation), name=str(score)) + "\n")

    config.dsRNA.predicted.parent.mkdir(parents=True, exist_ok=True)
    with open(config.dsRNA.predicted, 'wb') as stream:
        pickle.dump(pkl, stream, protocol=pickle.HIGHEST_PROTOCOL)

    with tempfile.NamedTemporaryFile() as tmp:
        with open(tmp.name, 'w') as stream:
            stream.writelines(track)
        utils.bed.tbindex(BedTool(tmp.name).sort(), config.dsRNA.predicted.with_suffix(".bed.gz"))
