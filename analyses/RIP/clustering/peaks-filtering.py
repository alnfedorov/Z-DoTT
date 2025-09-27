import pandas as pd
from pybedtools import BedTool, Interval

import ld

pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


def job(config: ld.Config):
    cache = config.peaks.filtered.with_suffix(".pkl")
    if not cache.exists():
        # Load all prefiltered peaks
        prefiltered = pd.read_pickle(config.peaks.prefiltered.with_suffix(".pkl"))

        # Annotated dsRNAs and replication
        prefiltered = ld.features.dsRNA(prefiltered, config)
        prefiltered = ld.features.replication(prefiltered, config, which='filtered')
        prefiltered.to_pickle(cache)

    prefiltered = pd.read_pickle(cache)

    mask = (prefiltered['dsRNA[Individual]'] == 0) & (prefiltered['Replication'] < config.peaks.strong_min_replication)
    prefiltered.loc[mask, 'auto-excluded'] = True
    filtered = prefiltered[~prefiltered['auto-excluded']].copy()

    filtered = ld.features.closest_neighbor_and_merged_length(filtered)
    filtered = filtered[(filtered['Merged length'] >= config.peaks.min_length) | (filtered['dsRNA[Individual]'] > 0)]

    # Select peaks passing the final filter and save them
    bed = [
        Interval(contig, start, end, strand=str(orientation))
        for contig, start, end, orientation in
        filtered[['contig', 'start', 'end', 'orientation']].itertuples(index=False, name=None)
    ]
    BedTool(bed).sort().saveas(config.peaks.filtered)


for config in ld.Config.load():
    job(config)
