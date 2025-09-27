import pickle
from collections import defaultdict

from biobit.core.loc import Strand, Interval
from biobit.toolkit.annotome.transcriptome import Location
from joblib import delayed, Parallel

import ld
import utils
from stories.annotation import TranscriptionBoundaries


def job(assembly: str):
    assembly = utils.assembly.get(name=assembly)

    print(f"Processing {assembly.name}...")

    boundaries = defaultdict(set)
    for seqid, size in assembly.seqid.sizes().items():
        for strand in Strand.Forward, Strand.Reverse:
            boundaries[seqid, strand].add((0, "seq-start"))
            boundaries[seqid, strand].add((size, "seq-end"))

    selected: list[tuple[Location, list[Interval]]] = []

    # GENCODE annotation
    gencode = ld.load.gencode(assembly.name)
    for rna in gencode.rnas.values():
        if ld.filters.is_transcription_boundary(rna):
            selected.append((rna.loc, rna.exons))

    # RefSeq annotation (same logic as above)
    for rna in assembly.refseq.load().rnas.values():
        if rna.ind.startswith("rna-NR_") or rna.ind.startswith("rna-NM_"):
            if ld.filters.is_transcription_boundary(rna):
                selected.append((rna.loc, rna.exons))

    # Process the selected annotations
    for loc, exons in selected:
        match loc.strand:
            case Strand.Forward:
                boundaries[loc.seqid, Strand.Forward].add((loc.start, "tss"))
                boundaries[loc.seqid, Strand.Forward].add((loc.end, "tes"))
                for prv, nxt in zip(exons[:-1], exons[1:]):
                    if prv.end != nxt.start:
                        boundaries[loc.seqid, Strand.Forward].add((prv.end, "donor"))
                        boundaries[loc.seqid, Strand.Forward].add((nxt.start, "acceptor"))
            case Strand.Reverse:
                boundaries[loc.seqid, Strand.Reverse].add((loc.end, "tss"))
                boundaries[loc.seqid, Strand.Reverse].add((loc.start, "tes"))
                for prv, nxt in zip(exons[:-1], exons[1:]):
                    if prv.end != nxt.start:
                        boundaries[loc.seqid, Strand.Reverse].add((prv.end, "acceptor"))
                        boundaries[loc.seqid, Strand.Reverse].add((nxt.start, "donor"))
            case _:
                raise ValueError(loc.strand)

    saveto = ld.paths.transcription_boundaries.pkl[assembly.name]
    saveto.parent.mkdir(parents=True, exist_ok=True)

    # Create the index
    index = {k: TranscriptionBoundaries(v) for k, v in boundaries.items()}
    with open(saveto, 'wb') as stream:
        pickle.dump(index, stream, protocol=pickle.HIGHEST_PROTOCOL)

    # Save as BED
    bed = utils.bed.blocks.from_iterable((
        (seqid, [Interval(pos, pos + 1)], strand.symbol(), name, None)
        for (seqid, strand), values in boundaries.items() for pos, name in values
    ))
    utils.bed.tbindex(bed, saveto.with_suffix(".bed.gz"))


# Save the result
Parallel(n_jobs=-1)(delayed(job)(assembly) for assembly in ("GRCm39", "CHM13v2"))
