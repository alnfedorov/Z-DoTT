from pathlib import Path
from typing import Any

from ..GRCm39.assembly import GRCm39
from ..HSV1.assembly import HSV1
from ..IAV.assembly import IAV
from ..resources import Assembly


class zdhs(Assembly):
    def __init__(
            self, fasta: Path, gff3: Path, hsv1: HSV1, iav: IAV, grcm39: GRCm39
    ):
        super().__init__()
        self.fasta = fasta
        self.gff3 = gff3
        self.hsv1 = hsv1
        self.iav = iav
        self.grcm39 = grcm39

    @property
    def name(self) -> str:
        return "zdmm"

    @property
    def organisms(self) -> frozenset[str]:
        return frozenset([*self.hsv1.organisms, *self.iav.organisms, *self.grcm39.organisms])

    def seqsizes(self) -> dict[str, int]:
        return self.hsv1.seqsizes() | self.iav.seqsizes() | self.grcm39.seqsizes()

    def seqlabels(self) -> dict[str, str]:
        return self.hsv1.seqlabels() | self.iav.seqlabels() | self.grcm39.seqlabels()

    def annotations(self) -> frozenset[str]:
        return frozenset([
            "REDIportal", "RepeatMasker", "RepeatMasker classification", "GFF3", "FASTA",
        ])

    def files(self, annotation: str) -> tuple[Path, ...]:
        match annotation:
            case "REDIportal" | "RepeatMasker" | "RepeatMasker classification":
                return self.grcm39.files(annotation)
            case "GFF3":
                return (self.gff3,)
            case "FASTA":
                return (self.fasta,)
            case _:
                raise ValueError(f"Unknown annotation: {annotation}")

    def load(self, annotation: str) -> Any:
        match annotation:
            case "REDIportal" | "RepeatMasker" | "RepeatMasker classification":
                return self.grcm39.load(annotation)
            case "GFF3" | "FASTA":
                raise NotImplementedError(f"Loading {annotation} is unsupported")
            case _:
                raise ValueError(f"Unknown annotation: {annotation}")
