from pathlib import Path

import pandas as pd
import pybedtools
from biobit.core.loc import Orientation, Interval
from joblib import Parallel, delayed

import ld
import utils
from assemblies import GRCm39, CHM13v2


def job(assembly):
    def load(sample: Path):
        df = ld.utils.load.reat(
            sample, ld.c.thresholds.min_coverage, ld.c.thresholds.min_edits, ld.c.thresholds.freqthr
        )
        return df[['contig', 'trstrand', 'pos']]

    # 1. Load all real samples
    dfs = Parallel(n_jobs=-1)(
        delayed(load)(sample)
        for sample in (ld.c.reat.candidates / assembly.name).glob("*/*.csv.gz")
        if not sample.name.endswith("[merged].csv.gz")
    )

    # 2. Collapse & count replication
    df = pd.concat(dfs)
    groups = {}
    for (seqid, strand), partition in df.groupby(['contig', 'trstrand']):
        groups[seqid, Orientation(strand)] = partition['pos'].value_counts().to_dict()

    # 3. Load splice sites index & REDI sites
    gencode = assembly.gencode.load()
    ssites = ld.utils.load.splice_sites(gencode)
    redi = ld.utils.load.sites(assembly.rediportal)

    # 4. Collect the final set of filtered editing sites
    def job(
            seqid: str, orient: Orientation, positions: dict[int, int],
            in_ssites: list[Interval], in_redi: list[Interval]
    ):
        print(f"[{assembly.name}] {seqid}[{orient}]: {len(positions)} candidates")
        replicated, unreplicated = [], []
        for pos, observed in positions.items():
            it = Interval(pos, pos + 1)
            if observed >= ld.c.thresholds.min_samples:
                replicated.append(it)
            else:
                unreplicated.append(it)

        # Rescue unreplicated sites recorded in REDI
        replicated = replicated + Interval.overlap(unreplicated, in_redi)
        assert sum(x.len() for x in replicated) == sum(x.len() for x in Interval.merge(replicated))

        # Filter out sites near splice sites
        before = sum(x.len() for x in replicated)
        replicated = Interval.subtract(replicated, in_ssites)
        after = sum(x.len() for x in replicated)

        return seqid, orient, replicated, before - after

    filtered = Parallel(n_jobs=-1)(
        delayed(job)(seqid, orient, positions, ssites.get((seqid, orient), []), redi.get((seqid, orient), []))
        for (seqid, orient), positions in groups.items()
    )

    passed, bed, ssignored = {}, [], 0
    for seqid, orient, replicated, ignored in filtered:
        if (seqid, orient) not in passed:
            passed[seqid, orient] = set()
        for site in replicated:
            assert site.len() == 1
            passed[seqid, orient].add(site.start)
        bed.extend([pybedtools.Interval(seqid, it.start, it.end, strand=orient.symbol()) for it in replicated])
        ssignored += ignored

    print(f"[{assembly.name}] Ignored sites near splice-sites: {ssignored}")

    # 5. Save all passed sites as BED file
    saveto = ld.c.tracks.all_passed / f"{assembly.name}.bed.gz"
    saveto.parent.mkdir(parents=True, exist_ok=True)

    utils.bed.tbindex(pybedtools.BedTool(bed).sort(), saveto)

    # 6. Filter REAT tables / Save individual tracks
    def filter(sample: Path, saveto: Path, trackto: Path):
        df = ld.utils.load.reat(
            sample, ld.c.thresholds.min_coverage, ld.c.thresholds.min_edits, ld.c.thresholds.freqthr
        )
        df['orientation'] = df['trstrand'].apply(Orientation)

        mask, intervals, remained = [], [], 0
        for contig, orient, pos, freq in zip(df['contig'], df['orientation'], df['pos'], df['freq']):
            if pos in passed.get((contig, orient), set()):
                mask.append(True)
                remained += 1
                intervals.append(pybedtools.Interval(
                    contig, pos, pos + 1, strand=orient.symbol(), score=int(freq * 1_000)
                ))
            else:
                mask.append(False)
        df[mask].to_csv(saveto, index=False)
        utils.bed.tbindex(pybedtools.BedTool(intervals).sort(), trackto)
        print(f"[{assembly.name}] {sample.with_suffix('').name} passed sites: {remained}({remained / len(df):.1%})")

    workload = []
    for dataset in (ld.c.reat.candidates / assembly.name).iterdir():
        saveto = ld.c.reat.filtered / assembly.name / dataset.name
        saveto.mkdir(parents=True, exist_ok=True)

        trackto = ld.c.tracks.expwise / assembly.name / dataset.name
        trackto.mkdir(parents=True, exist_ok=True)

        for sample in dataset.glob("*.csv.gz"):
            workload.append((
                sample,
                saveto / sample.name,
                trackto / (sample.name.split('.')[0] + ".bed.gz")
            ))

    Parallel(n_jobs=-1)(delayed(filter)(*args) for args in workload)


Parallel(n_jobs=-1, backend='sequential')(delayed(job)(assembly) for assembly in [GRCm39, CHM13v2])
