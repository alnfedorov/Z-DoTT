import pickle

import pandas as pd
from biobit.core.loc import Strand
from biobit.toolkit.libnorm import MedianOfRatiosNormalization
from pybedtools.cbedtools import defaultdict

from .ld import FRAGMENTS, SERIES, BINS, VIRAL_RNA_LOAD, NormBin


def median_of_ratios() -> dict[str, MedianOfRatiosNormalization]:
    # Load all scaling bins
    with open(FRAGMENTS, 'rb') as stream:
        _scbins = pickle.load(stream)

    scaling = {}
    for assembly, bins in _scbins.items():
        bins = bins.set_index('source').T
        scaling[assembly] = MedianOfRatiosNormalization(bins, minval=10)
    return scaling


def load_counts() -> dict[str, pd.DataFrame]:
    with open(FRAGMENTS, 'rb') as stream:
        counts = pickle.load(stream)
    return counts


def load_vRNA_load() -> pd.DataFrame:
    with open(VIRAL_RNA_LOAD, 'rb') as stream:
        vRNA_LOAD = pickle.load(stream)
    return vRNA_LOAD


def load_bins() -> dict[str, dict[str, NormBin]]:
    with open(BINS, 'rb') as stream:
        bins: dict[str, dict[tuple[str, Strand], tuple[NormBin, ...]]] = pickle.load(stream)

    flatten = defaultdict(list)
    for assembly, partitions in bins.items():
        for bins in partitions.values():
            flatten[assembly].extend(bins)

    flatten = {
        assembly: {bin.ind: bin for bin in bins}
        for assembly, bins in flatten.items()
    }
    return flatten
