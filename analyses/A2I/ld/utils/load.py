from collections import defaultdict
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from biobit.core.loc import Orientation, Interval
from biobit.toolkit import seqproj
from biobit.toolkit.annotome import Annotome
from joblib import Parallel, delayed
from pybedtools import BedTool


def splice_sites(annotation: Annotome, offset: int = 15) -> dict[tuple[str, Orientation], list[Interval]]:
    ss = defaultdict(list)

    for rna in annotation.rnas.values():
        for prv, nxt in zip(rna.exons[:-1], rna.exons[1:]):
            if rna.loc.strand == "+":
                ss[rna.loc.seqid, Orientation.Forward].append(Interval(prv.end, prv.end + offset))
                ss[rna.loc.seqid, Orientation.Forward].append(Interval(nxt.start - offset, nxt.start))
            else:
                assert rna.loc.strand == "-"
                ss[rna.loc.seqid, Orientation.Reverse].append(Interval(prv.end, prv.end + offset))
                ss[rna.loc.seqid, Orientation.Reverse].append(Interval(nxt.start - offset, nxt.start))

    # Merge overlapping intervals
    ss = {key: Interval.merge(intervals) for key, intervals in ss.items()}
    return ss


def sites(path: Path) -> dict[tuple[str, Orientation], list[Interval]]:
    index = defaultdict(list)
    for site in BedTool(path):
        index[site.chrom, Orientation(site.strand)].append(Interval(site.start, site.end))
    index = {key: Interval.merge(intervals) for key, intervals in index.items()}
    return index


def reat(sample: Path, min_coverage: int, min_edits: int, freqthr: tuple[float, float]):
    df = pd.read_csv(sample, dtype={
        "contig": str, "trstrand": str, "pos": np.uint32, "refnuc": np.uint32, "misnuc": np.uint32
    })
    df['freq'] = df['misnuc'] / (df['refnuc'] + df['misnuc'])

    mask = (df['freq'] >= freqthr[0]) & (df['freq'] <= freqthr[1])
    mask &= (df['refnuc'] + df['misnuc']) >= min_coverage
    mask &= df['misnuc'] >= min_edits

    return df[mask]


def filtered(sample: Path):
    return pd.read_csv(sample, dtype={
        "contig": str, "trstrand": str, "pos": np.uint32, "refnuc": np.uint32, "misnuc": np.uint32,
        "freq": np.float32
    })


def annotated(
        root: Path, experiments: dict[str, list[seqproj.Experiment]],
        metainfo: Callable[[seqproj.Experiment], dict[str, str]],
        locmapping: dict[str, str] | None = None, pivot: list[str] = None
) -> pd.DataFrame:
    def loadit(prjind: str, expind: str, info: dict[str, str]):
        path = root / prjind / f"{expind}.csv.gz"
        df = pd.read_csv(path, dtype={
            "contig": str, "trstrand": str, "pos": np.uint32, "refnuc": np.uint32, "misnuc": np.uint32,
            "freq": np.float32, "location": str
        })
        if locmapping:
            mapped = df['location'].apply(lambda x: locmapping.get(x, pd.NA))
            skipped = mapped.isnull()
            if skipped.any():
                print(f"Skipped locations: {df.loc[skipped, 'location'].value_counts().to_dict()}")
            df['location'] = mapped
            df = df[~skipped].copy()

        for key, value in info.items():
            df[key] = value
        return df

    workload = []
    for project, experiments in experiments.items():
        prjind = project.replace(" ", "_").replace("/", "_")
        for exp in experiments:
            expind = exp.ind.replace(" ", "_").replace("/", "_")
            workload.append((prjind, expind, metainfo(exp)))

    results = Parallel(n_jobs=-1)(delayed(loadit)(prjind, expind, meta) for prjind, expind, meta in workload)
    df = pd.concat(results)

    if pivot is not None:
        df['total'] = 1
        total = df.groupby(pivot)['total'].sum().to_dict()
        df['sites'] = [w / total[key] for w, key in zip(df['total'], zip(*[df[col] for col in pivot]))]

        # Collapse locations & normalize
        groups = [*pivot, 'location']
        df = df.groupby(groups)[['sites', 'total']].sum().reset_index()
        df['sites'] *= 100
    return df
