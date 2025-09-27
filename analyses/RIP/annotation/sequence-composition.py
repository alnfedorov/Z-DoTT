import pickle
from collections import defaultdict
from pathlib import Path
from typing import Any

import intervaltree
import pandas as pd
from biobit.core.loc import Orientation
from biobit.deprecated.gindex import GenomicIndex
from joblib import Parallel, delayed

import ld
import utils
from assemblies import GRCm39, CHM13v2
from stories.RIP import annotation, pcalling


def index_repeats(name: str, repmasker: Path) -> tuple[str, GenomicIndex]:
    itrees = defaultdict(lambda: intervaltree.IntervalTree())
    bed = pd.read_csv(
        repmasker, sep='\t', header=None, names=['contig', 'start', 'end', 'name', 'score', 'strand'],
        dtype={'contig': str, 'start': int, 'end': int, 'name': str, 'score': str, 'strand': str}
    )
    for contig, start, end, repname, _, strand in bed.itertuples(index=False):
        itrees[(contig, Orientation.Dual)].addi(start, end, data=(repname, strand))
    return name, GenomicIndex(itrees)


# Index the repeat annotations
_index = Parallel(n_jobs=-1, verbose=100)(
    delayed(index_repeats)(assembly.name, assembly.repmasker) for assembly in [GRCm39, CHM13v2]
)
INDEX: dict[str, GenomicIndex] = dict(_index)

# Load the signal cache
SIGNAL = ld.cache.signal.load()


def resolve_sequence(repmasker, scores: dict[Any, float]) -> dict[str, float]:
    if not scores:
        return {"Undefined": 1}

    minval = min(scores.values()) - 1e-6
    scores = {k: v - minval for k, v in scores.items()}

    fractions, total = defaultdict(float), sum(scores.values())
    assert total > 0, scores
    for (left, right), v in scores.items():
        fractions[left, right] += v / total
    assert abs(sum(fractions.values()) - 1) < 1e-3, (fractions, scores)
    return fractions


def job(config: annotation.Config, cmp: pcalling.Config):
    host = utils.assembly.get(organism=cmp.host)
    index = INDEX[host.name]
    signal = SIGNAL[config.ind]

    results = {}
    for partition in config.elements:
        cache = signal[partition.contig, partition.orientation]
        scores = defaultdict(int)

        # Calculate scores for individual peaks
        for peak in partition.peaks:
            score = cache[cmp.project, cmp.ind][peak]
            for segment, annotation in index.overlap(partition.contig, Orientation.Dual, rng=peak).to_steps():
                if not annotation:
                    annotation = {("Repeat-free", ".")}
                for anno in annotation:
                    # scores[anno] += segment.len() / peak.len() * score
                    scores[anno] = max(scores[anno], score)

        for invrep in partition.invrep:
            for arm in invrep.segments:
                armlength = arm.left.len()
                score = cache[cmp.project, cmp.ind][arm.left] + cache[cmp.project, cmp.ind][arm.right]

                # Annotate left/right InvRepeat arms
                left = index.overlap(partition.contig, Orientation.Dual, rng=arm.left).to_steps()
                left = list(zip(left.boundaries, left.annotations))

                # And flip the right arm
                right = index.overlap(partition.contig, Orientation.Dual, rng=arm.right).to_steps()
                right = list(zip(right.boundaries[::-1], right.annotations[::-1]))

                lind, rind, covered = 0, 0, 0
                while True:
                    # Calculate current offsets in the segment coordinates
                    loff = (left[lind][0].start - arm.left.start, left[lind][0].end - arm.left.start)
                    roff = (arm.right.end - right[rind][0].end, arm.right.end - right[rind][0].start)

                    # assert loff[0] == roff[0]
                    offset = (max(loff[0], roff[0]), min(loff[1], roff[1]))
                    offlength = offset[1] - offset[0]

                    # Add the score to each annotation combination
                    lanno = left[lind][1] if left[lind][1] else {("Repeat-free", ".")}
                    ranno = right[rind][1] if right[rind][1] else {("Repeat-free", ".")}
                    for l in lanno:
                        for r in ranno:
                            # scores[l, r] += score * offlength / armlength
                            scores[l, r] = max(scores[l, r], score)

                    covered += offset[1] - offset[0]
                    # Move pointers if needed
                    if loff[1] == offset[1]:
                        lind += 1
                    if roff[1] == offset[1]:
                        rind += 1

                    if lind == len(left) and rind == len(right):
                        break
                assert covered == armlength

        results[partition.ind] = resolve_sequence(host.repcls, scores)
    return config.ind, cmp.ind, results


configs = annotation.Config.load()
_results = Parallel(n_jobs=-1, verbose=100, backend='multiprocessing')(
    delayed(job)(config, cmp) for config in configs for cmp in config.comparisons
)

# Group the results and save to the cache
results = defaultdict(dict)
for ind, cmp, res in _results:
    results[ind][cmp] = res

for conf in configs:
    conf.sequences.parent.mkdir(parents=True, exist_ok=True)
    with open(conf.sequences, 'wb') as stream:
        pickle.dump(results[conf.ind], stream)
