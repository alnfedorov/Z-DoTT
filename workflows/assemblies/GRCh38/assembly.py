from pathlib import Path
from typing import Any

from . import seqid
from ..resources import Assembly


class GRCh38(Assembly):
    def __init__(self, fasta: Path, gencode_gff3: Path):
        super().__init__()
        self.fasta = fasta
        self.gencode_gff3 = gencode_gff3

    @property
    def name(self) -> str:
        return "GRCh38"

    @property
    def organisms(self) -> frozenset[str]:
        return frozenset(["Homo sapiens"])

    def seqsizes(self) -> dict[str, int]:
        return seqid.sizes()

    def annotations(self) -> frozenset[str]:
        return frozenset(["GENCODE", "FASTA"])

    def files(self, annotation: str) -> tuple[Path, ...]:
        match annotation:
            case "GENCODE":
                return (self.gencode_gff3,)
            case "FASTA":
                return (self.fasta,)
            case _:
                raise ValueError(f"Unknown annotation: {annotation}")

    def load(self, annotation: str) -> Any:
        raise NotImplementedError("Loading of annotations is not implemented for GRCh38.")
