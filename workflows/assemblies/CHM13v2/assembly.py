from pathlib import Path
from typing import Any

from pybedtools import BedTool

from . import seqid
from ..resources import Assembly


class CHM13v2(Assembly):
    def __init__(
            self, fasta: Path, gencode_liftoff_gff3: Path, refseq_gff3: Path, rediportal: Path,
            repmasker: Path, repmasker_classification: Path
    ):
        super().__init__()
        self.fasta = fasta
        self.gencode_liftoff_gff3 = gencode_liftoff_gff3
        self.refseq_gff3 = refseq_gff3
        self.rediportal = rediportal
        self.repmasker = repmasker
        self.repmasker_classification = repmasker_classification

    @property
    def name(self) -> str:
        return "CHM13v2"

    @property
    def organisms(self) -> frozenset[str]:
        return frozenset(["Homo sapiens"])

    def seqsizes(self) -> dict[str, int]:
        return seqid.sizes()

    def annotations(self) -> frozenset[str]:
        return frozenset([
            "REDIportal", "RepeatMasker", "RepeatMasker classification", "GENCODE-Liftoff", "RefSeq", "FASTA",
        ])

    def files(self, annotation: str) -> tuple[Path, ...]:
        match annotation:
            case "REDIportal":
                return (self.rediportal,)
            case "RepeatMasker":
                return (self.repmasker,)
            case "RepeatMasker classification":
                return (self.repmasker_classification,)
            case "GENCODE-Liftoff":
                return (self.gencode_liftoff_gff3,)
            case "RefSeq":
                return (self.refseq_gff3,)
            case "FASTA":
                return (self.fasta,)
            case _:
                raise ValueError(f"Unknown annotation: {annotation}")

    def load(self, annotation: str) -> Any:
        match annotation:
            case "REDIportal":
                return BedTool(self.rediportal)
            case "RepeatMasker":
                return BedTool(self.repmasker)
            case "RepeatMasker classification":
                return self.repmasker_classification
            case "GENCODE-Liftoff" | "RefSeq" | "FASTA":
                raise NotImplementedError(f"Loading {annotation} is unsupported")
            case _:
                raise ValueError(f"Unknown annotation: {annotation}")
