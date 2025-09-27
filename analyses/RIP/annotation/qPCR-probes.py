import pickle
from collections import defaultdict

import intervaltree as it
from biobit.core.loc import Orientation, Interval
from joblib import Parallel, delayed
from tqdm import tqdm

from stories.RIP import annotation

# Rough coordinates of the qPCR probes in each set of experiments
PROBES = {
    "Z-RNA[CHM13v2]": {
        "C6orf62": ("chr6", "-", Interval(24539443, 24547846)),
        "RUVBL2": ("chr19", "+", Interval(52005293, 52030293)),
        "Inverted ZNFs": ("chr19", "-", Interval(61511932, 61531665)),
        "RABGGTB": ("chr1", "+", Interval(75659652, 75672126)),
        "ROCK1": ("chr18", "-", Interval(21119286, 21124548)),
        "NCL": ("chr2", "-", Interval(231932148, 231934354)),
        "KBTBD2": ("chr7", "-", Interval(32970998, 32986006)),
    },
    "Z-RNA[GRCm39]": {
        "Btbd3": ("chr2", "+", Interval(138174403, 138192278)),
        "Ptbp1": ("chr10", "+", Interval(79698893, 79712686)),
        "Gtpbp4": ("chr13", "-", Interval(8988283, 9014023)),
        "H2ac18/19(+)": ("chr3", "+", Interval(96147973, 96152604)),
        "H2ac18/19(-)": ("chr3", "-", Interval(96147973, 96152604)),
        "Nabp1": ("chr1", "-", Interval(51470634, 51483545)),
        "Hmga1": ("chr17", "-", Interval(27743519, 27752636)),
        "U1 locus": ("chr3", "-", Interval(96394473, 96396018)),
        "Rock1": ("chr18", "-", Interval(10040233, 10059489)),
        "Vbp1": ("chrX", "+", Interval(74600405, 74603219)),
        "Gm14419": ("chr2", "+", Interval(176725754, 176730032)),
        "Fam24a/b": ("chr7", "-", Interval(130927631, 130929229)),
        "Gnpda1": ("chr18", "-", Interval(38457830, 38458499)),
        "Anapc4": ("chr5", "-", Interval(52987842, 52990626)),
        "Pcna": ("chr2", "-", Interval(132084025, 132085391)),
        "Haus2": ("chr2", "+", Interval(120450019, 120454661)),
        "Tlcd1": ("chr11", "+", Interval(78074357, 78078886)),
    }
}


def job(config: annotation.Config):
    # Make the probes index
    probes = defaultdict(it.IntervalTree)
    for name, (seqid, orientation, interval) in PROBES[config.ind].items():
        probes[seqid, Orientation(orientation)].addi(interval.start, interval.end, data=name)

    results = {}
    for partition in tqdm(config.elements):
        index = probes[partition.contig, partition.orientation]

        overlap = set()
        for peak in partition.peaks:
            overlap |= index.overlap(peak.start, peak.end)
        for invrep in partition.invrep:
            for seqrng in invrep.seqranges():
                overlap |= index.overlap(seqrng.start, seqrng.end)

        if not overlap:
            results[partition.ind] = None
        elif len(overlap) == 1:
            results[partition.ind] = overlap.pop().data
        else:
            raise ValueError(f"Multiple probes overlap the same region: {overlap}")

    expected = set(PROBES[config.ind])
    observed = set(results.values()) - {None}
    if expected != observed:
        raise ValueError(f"Expected probes {expected} but got {observed}. Missing: {expected - observed}")

    config.probes.parent.mkdir(parents=True, exist_ok=True)
    with open(config.probes, 'wb') as stream:
        pickle.dump(results, stream)


Parallel(n_jobs=-1, verbose=100)(delayed(job)(config) for config in annotation.Config.load())
