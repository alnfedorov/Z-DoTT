import pickle

from biobit import io
from biobit.core.loc import Strand
from biobit.core.ngs import Layout
from biobit.toolkit import countit, nfcore

import ld
from assemblies import CHM13v2, HSV1, IAV, GRCm39
from stories.normalization.ld import NormBin


def reader(path: str, layout: Layout) -> io.bam.Reader:
    assert isinstance(layout, Layout.Paired), f"Unsupported layout: {layout} ({path})"
    return io.bam.Reader(path, inflags=3, exflags=2572, minmapq=0)


with open(ld.BINS, 'rb') as stream:
    BINS: dict[str, dict[tuple[str, Strand], tuple[NormBin, ...]]] = pickle.load(stream)


# Setup a CountIt engine for each assembly
def setup(assembly) -> countit.rigid.Engine[str]:
    bins = BINS[assembly.name]

    partitions, annotations = set(), []
    for (seqid, strand), bins in bins.items():
        partitions.add(seqid)
        for bin in bins:
            annotations.append((bin.ind, [(seqid, strand.to_orientation(), list(bin.intervals))]))

    sizes = assembly.seqid.sizes() | IAV.segments | HSV1.segments
    partitions = [(contig, (0, size)) for contig, size in sizes.items()]

    engine = countit.rigid.Engine.builder().set_threads(-1).add_elements(annotations).add_partitions(partitions).build()
    return engine


allcounts, allmetrics = {}, {}
for assembly, datasets in ld.SERIES.items():
    if not datasets:
        continue

    assembly = {"GRCm39": GRCm39, "CHM13v2": CHM13v2}[assembly]
    engine = setup(assembly)

    # Run the counting
    sources = []
    for prj in datasets:
        for exp_ind, (bam, layout) in nfcore.rnaseq.extract.bams(prj, factory=reader).items():
            sources.append(((prj.ind, exp_ind), bam, layout))
    resolution = countit.rigid.resolution.OverlapWeighted()
    result = engine.run(sources, resolution)

    cnts, mts = countit.utils.result_to_pandas(result)
    mts["segment"] = mts["segment"].apply(str)

    assert assembly.name not in allcounts and assembly.name not in allmetrics
    allcounts[assembly.name] = cnts
    allmetrics[assembly.name] = mts

for element, saveto in (allcounts, ld.FRAGMENTS), (allmetrics, ld.METRICS):
    saveto.parent.mkdir(parents=True, exist_ok=True)
    with open(saveto, "wb") as f:
        pickle.dump(element, f)
