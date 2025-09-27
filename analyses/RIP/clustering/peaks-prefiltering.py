from multiprocessing import Pool, cpu_count

import pandas as pd
from pybedtools import BedTool, Interval

import ld

pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


def job(config: ld.Config):
    cache = config.peaks.universe.with_suffix(".pkl")
    if not cache.exists():
        # Derive all covered regions to use as natural insulators for the repeto groups
        ld.features.derive_covered_regions(config)

        # Create the peaks universe: all observed peak pieces
        universe = ld.features.build_universe(config)

        # Derive all the relevant annotations
        universe = ld.features.replication(universe, config, which='filtered')
        universe = ld.features.closest_neighbor_and_merged_length(universe)

        pool = Pool(cpu_count())
        for ind, handle in enumerate([
            pool.apply_async(ld.features.repeats, args=(universe, config)),
            pool.apply_async(ld.features.genomic_regions, args=(universe, config)),
            pool.apply_async(ld.features.curated_regions, args=(universe, config)),
            pool.apply_async(ld.features.editing_sites, args=(universe, config)),
        ]):
            handle.wait()
            assert handle.successful(), f"Failing on {ind}"
            for col, values in handle.get().items():
                universe[col] = values

        cache.parent.mkdir(parents=True, exist_ok=True)
        with open(cache, 'wb') as stream:
            universe.to_pickle(stream)
    else:
        with open(cache, 'rb') as stream:
            universe = pd.read_pickle(stream)

    # Annotated curated regions
    impossible = universe['Curated[included]'] & universe['Curated[excluded]']
    assert impossible.sum() == 0, universe[impossible]

    # All regions are included by default
    universe['auto-excluded'] = False

    # Ignore viral contigs
    universe.loc[universe['Region'].isna() | (universe['Region'] == 'Not host'), 'auto-excluded'] = True

    # rRNA-like sequences
    mask = (universe['RepeatMasker'] == 'SSU-rRNA_Hsa') | (universe['RepeatMasker'] == 'LSU-rRNA_Hsa')
    universe.loc[mask, 'auto-excluded'] = True

    # Homo polymers
    for x in "A", "T", "G", "C":
        universe.loc[universe['RepeatMasker'] == f'({x})n', 'auto-excluded'] = True

    # Low replication
    universe['auto-excluded'] |= universe['Replication'] < config.peaks.relaxed_min_replication

    # Very low length
    universe.loc[universe['Merged length'] <= config.peaks.min_length, 'auto-excluded'] = True

    # Keep edited regions
    include = (universe['Editing sites'] >= config.peaks.relaxed_min_editing_sites)  # | universe['Curated[included]']
    universe.loc[include, 'auto-excluded'] = False

    # Curated regions
    universe['auto-excluded'] |= universe['Curated[excluded]']

    # Select peaks passing the initial filter and save them
    passed = universe[~universe['auto-excluded']]
    for peaks, saveto in (universe, config.peaks.universe), (passed, config.peaks.prefiltered):
        bed = [
            Interval(contig, start, end, strand=str(orientation), name=f"N={replication}")
            for contig, start, end, orientation, replication in
            peaks[['contig', 'start', 'end', 'orientation', 'Replication']].itertuples(index=False, name=None)
        ]
        BedTool(bed).sort().saveas(saveto)

    # Save prefiltered meta for the later use
    passed.to_pickle(config.peaks.prefiltered.with_suffix(".pkl"))


for config in ld.Config.load():
    job(config)
