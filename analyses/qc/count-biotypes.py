import pickle
from collections import defaultdict

import pandas as pd
from biobit import io
from biobit.core.ngs import Layout
from biobit.toolkit import countit, nfcore

import ld
from assemblies import GRCm39, CHM13v2, HSV1, IAV
from stories import annotation


def reader(path: str, layout: Layout) -> io.bam.Reader:
    assert isinstance(layout, Layout.Paired), f"Unsupported layout: {layout} ({path})"
    return io.bam.Reader(path, inflags=3, exflags=2572, minmapq=0)


def setup(assembly) -> countit.rigid.Engine[str]:
    ant = annotation.load.resolved_annotation(assembly.name)

    elements = [
        (key, [(contig, orient, segments) for (contig, orient), segments in data.items()])
        for key, data in ant.items()
    ]

    seqids = assembly.seqid.sizes() | HSV1.segments | IAV.segments
    partitions = [(contig, (0, size)) for contig, size in seqids.items()]

    # Create the CountIt engine
    engine = countit.rigid.Engine.builder().set_threads(-1).add_elements(elements).add_partitions(partitions).build()

    return engine


counts, metrics = [], []
for assembly, datasets in ld.BIOTYPES.series.items():
    if not datasets:
        continue

    assembly = {"GRCm39": GRCm39, "CHM13v2": CHM13v2}[assembly]
    ctit = setup(assembly)

    # Run the counting
    sources = []
    for prj in datasets:
        for exp_ind, (bam, layout) in nfcore.rnaseq.extract.bams(prj, factory=reader).items():
            sources.append(((prj.ind, exp_ind), bam, layout))
    resolution = countit.rigid.resolution.OverlapWeighted()

    result = ctit.run(sources, resolution)
    cnts, mts = countit.utils.result_to_pandas(result)

    cnts["Assembly"], mts["Assembly"] = assembly.name, assembly.name
    counts.append(cnts)
    metrics.append(mts)

ld.BIOTYPES.METRICS.parent.mkdir(parents=True, exist_ok=True)
allmetrics = pd.concat(metrics)
allmetrics["segment"] = allmetrics["segment"].apply(str)
allmetrics.to_pickle(ld.BIOTYPES.METRICS)

ld.BIOTYPES.RAW_FRAGMENTS.parent.mkdir(parents=True, exist_ok=True)
with open(ld.BIOTYPES.RAW_FRAGMENTS, "wb") as f:
    pickle.dump(counts, f)

collapsed = []
for cnts in counts:
    groups = defaultdict(list)
    for col in cnts.select_dtypes("number").columns:
        groups[col[0]].append(col)
    for name, cols in groups.items():
        cnts[name] = cnts[cols].sum(axis=1)
        cnts = cnts.drop(columns=cols)
    collapsed.append(cnts)

collapsed = pd.concat(collapsed).fillna(0)

ld.BIOTYPES.COLLAPSED_FRAGMENTS.parent.mkdir(parents=True, exist_ok=True)
collapsed.to_pickle(ld.BIOTYPES.COLLAPSED_FRAGMENTS)
