import pickle
from collections import defaultdict

from biobit.core.loc import Strand, ChainInterval, Interval
from joblib import Parallel, delayed
from pybedtools import BedTool

import ld
import utils
from stories import annotation
from stories.annotation import RNACore, TranscriptionBoundaries


def job(
        assembly: str, seqid: str, strand: Strand, rna_cores: list[RNACore], boundaries: TranscriptionBoundaries
):
    chains = []
    allcovered = []

    # Canonical RNAs extended at 3' and 5' ends up to the longest annotated TES/TSS
    upstream, downstream = ('left', 'right') if strand == Strand.Forward else ('right', 'left')
    for rna in rna_cores:
        exons = rna.exons
        start, end = exons[0].start, exons[-1].end

        # Extend the RNA upstream until the next non-TSS region
        window, key = boundaries.window(start if strand == Strand.Forward else end, upstream)
        while key == 'tss':
            window, key = boundaries.window(window.start if strand == Strand.Forward else window.end, upstream)
        tss = window.end if strand == Strand.Forward else window.start

        # Extend the RNA downstream until the next non-TES region
        window, key = boundaries.window(end if strand == Strand.Forward else start, downstream)
        while key == 'tes':
            window, key = boundaries.window(window.end if strand == Strand.Forward else window.start, downstream)
        tes = window.start if strand == Strand.Forward else window.end

        newstart, newend = (tss, tes) if strand == Strand.Forward else (tes, tss)
        assert newstart <= start < end <= newend, ([newstart, newend], [start, end])
        exons = Interval.merge([Interval(newstart, exons[0].end), *exons[1:-1], Interval(exons[-1].start, newend)])

        chains.append(ChainInterval(exons))
        allcovered.extend(chains[-1])

    # Add each yet uncovered region as a separate chain
    allcovered = Interval.merge(allcovered)

    position, ind = 0, 0
    while ind < len(allcovered):
        assert position <= allcovered[ind].start
        if position == allcovered[ind].start:
            position = allcovered[ind].end
            ind += 1
            continue

        # Insert all transcriptionally independent regions between the current and the covered region
        nxt = boundaries.closest(position, 'right')[0]
        nxt = min(allcovered[ind].start, nxt)

        chains.append(ChainInterval([Interval(position, nxt)]))
        position = nxt

    while position != boundaries.boundaries[-1]:
        nxt = boundaries.closest(position, 'right')[0]
        chains.append(ChainInterval([Interval(position, nxt)]))
        position = nxt
    return assembly, seqid, strand, chains


rna_cores, boundaries = {}, {}
for assembly in "CHM13v2", "GRCm39":
    rna_cores[assembly] = annotation.load.rna_cores(assembly)
    boundaries[assembly] = annotation.load.transcription_boundaries(assembly)

result = Parallel(n_jobs=-1)(
    delayed(job)(
        assembly, seqid, strand, rna_cores[assembly].get((seqid, strand), []), boundaries[assembly][seqid, strand]
    )
    for assembly in boundaries
    for seqid, strand in boundaries[assembly]
)

parsed = defaultdict(dict)
bed = defaultdict(list)
for assembly, seqid, strand, chains in result:
    assert (seqid, strand) not in parsed[assembly]
    parsed[assembly][seqid, strand] = sorted(chains, key=lambda x: list(x)[0].start)

    for ind, chain in enumerate(chains):
        bed[assembly].append(utils.bed.blocks.make(
            seqid, chain, strand.to_orientation(), name=f"{ind}[{seqid}, {strand}]"
        ))

ld.RESULTS.mkdir(parents=True, exist_ok=True)
for assembly in bed:
    utils.bed.tbindex(
        BedTool(bed[assembly]).sort(),
        ld.RESULTS / f"{assembly}_models.bed.gz"
    )

ld.RNA_MODELS.parent.mkdir(parents=True, exist_ok=True)
with open(ld.RNA_MODELS, 'wb') as stream:
    pickle.dump(dict(parsed), stream, protocol=pickle.HIGHEST_PROTOCOL)
