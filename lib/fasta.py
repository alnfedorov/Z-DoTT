from functools import lru_cache
from pathlib import Path

import pysam


@lru_cache(maxsize=10)
def _fasta(fasta: Path) -> pysam.FastaFile:
    return pysam.FastaFile(fasta.as_posix())


@lru_cache(maxsize=10)
def contigs(fasta: Path) -> dict[str, int]:
    f = _fasta(fasta)
    return dict(zip(f.references, f.lengths))


def sequence(fasta: Path, contig: str, start: int, end: int, strand: str = "+") -> str:
    seq = _fasta(fasta).fetch(contig, start, end).upper()

    if strand == '-':
        seq = seq[::-1]
        seq = seq.translate(str.maketrans("ACGT", "TGCA"))
    return str(seq)
