import pickle
from typing import Any

import pandas as pd
import pybedtools
from biobit import io
from biobit.core.ngs import Layout
from biobit.toolkit import countit, nfcore

import ld
import utils
from stories.annotation import RNACore


def reader(path: str, layout: Layout) -> io.bam.Reader:
    assert isinstance(layout, Layout.Paired), f"Unsupported layout: {layout} ({path})"
    return io.bam.Reader(path, inflags=3, exflags=2572, minmapq=0)


with open(ld.RNA, 'rb') as stream:
    RNA = pickle.load(stream)


# Setup a CountIt engine for each assembly
def setup(assembly) -> tuple[countit.rigid.Engine[tuple[str, str]], dict[Any, int]]:
    annotation = []
    for seqid, rnas in RNA[assembly.name].items():
        for rna in rnas:  # type: RNACore
            orient = rna.strand.to_orientation()

            # Truncate aberrant regions
            rna = rna.truncate(
                max_read_through=ld.thr.scores.winsize.read_through,
                max_read_in=ld.thr.scores.winsize.read_in,
                max_divergent=ld.thr.scores.winsize.divergent
            )
            targets = []
            allexons = list(rna.exons)

            # Read-through transcription
            if rna.read_through.len() == ld.thr.scores.winsize.read_through:
                targets.append([("read-through", "reference", rna.gid), orient, allexons])
                targets.append([("read-through", "aberrant", rna.gid), orient, [rna.read_through]])

            # Read-in transcription
            if rna.read_in.len() == ld.thr.scores.winsize.read_in:
                targets.append([("read-in", "reference", rna.gid), orient, allexons])
                targets.append([("read-in", "aberrant", rna.gid), orient, [rna.read_in]])

            # Divergent transcription
            if rna.divergent.len() == ld.thr.scores.winsize.divergent:
                dorient = rna.strand.flipped().to_orientation()
                targets.append([("divergent", "reference", rna.gid), orient, allexons])
                targets.append([("divergent", "aberrant", rna.gid), dorient, [rna.divergent]])

            # Intronic transcription
            if rna.introns:
                for ind, intron in enumerate(rna.introns):
                    gid = f"{rna.gid}-{ind}"
                    targets.append([("intronic", "reference", gid), orient, [intron.donor]])
                    targets.append([("intronic", "aberrant", gid), orient, list(intron.intron)])

            for key, orientation, regions in targets:
                annotation.append((key, [(seqid, orientation, regions)]))

    partitions = [(seqid, (0, size)) for seqid, size in assembly.seqid.sizes().items()]

    # Save as a BED file
    bed = pybedtools.BedTool([
        utils.bed.blocks.make(contig, segments, str(orient), name="-".join(map(str, key)))
        for key, items in annotation
        for contig, orient, segments in items
    ])
    ld.counts.root.mkdir(parents=True, exist_ok=True)
    utils.bed.tbindex(bed, ld.counts.root / f"{assembly.name}.counting.bed.gz")

    # Extract the length of each region
    lengths = {}
    for key, items in annotation:
        length = 0
        for _, _, regions in items:
            length += sum(x.len() for x in regions)
        if key in lengths:
            assert lengths[key] == length, f"{key} has different lengths"
        lengths[key] = length

    # Create the CountIt engine
    engine = countit.rigid.Engine.builder() \
        .set_threads(-1) \
        .add_elements(annotation) \
        .add_partitions(partitions) \
        .build()
    return engine, lengths


counts, metrics, lengths = [], [], {}
for assembly, datasets in ld.comparisons.series.items():
    if not datasets:
        continue

    assembly = utils.assembly.get(name=assembly)
    ctit, lens = setup(assembly)

    # Record the lengths
    assert len(set(lengths) & set(lens)) == 0, "Duplicate length keys"
    lengths.update(lens)

    # Run the counting
    sources = []
    for prj in datasets:
        for exp_ind, (bam, layout) in nfcore.rnaseq.extract.bams(prj, factory=reader).items():
            sources.append(((prj.ind, exp_ind), bam, layout))
    resolution = countit.rigid.resolution.AnyOverlap()

    result = ctit.run(sources, resolution)
    cnts, mts = countit.utils.result_to_pandas(result)

    counts.append(cnts)
    metrics.append(mts)

# Save metrics
allmetrics = pd.concat(metrics)
with open(ld.counts.metrics, "wb") as f:
    pickle.dump(allmetrics, f)

# Postprocess counts
allcnts = pd.concat(counts).fillna(0)
allcnts = allcnts.set_index('source').T.reset_index()

allcnts['Category'] = allcnts['index'].apply(lambda x: x[0])
allcnts['Region'] = allcnts['index'].apply(lambda x: x[1])
allcnts['Ensembl ID'] = allcnts['index'].apply(lambda x: x[2])
allcnts['Length'] = allcnts['index'].apply(lambda x: lengths[x])

allcnts = allcnts.drop(columns=['index'])

# Save relevant regions
for region, saveto in [
    ("read-through", ld.counts.read_through),
    ("read-in", ld.counts.read_in),
    ("intronic", ld.counts.intronic),
    ("divergent", ld.counts.divergent)
]:
    df = allcnts[allcnts['Category'] == region].drop(columns='Category').copy()

    # Drop empty records
    mask = df.select_dtypes('number').sum(axis=1) > 1
    df = df[mask].groupby('Ensembl ID').filter(lambda x: len(x) == 2)

    df.to_pickle(saveto)
