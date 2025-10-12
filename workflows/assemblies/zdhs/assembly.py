from pathlib import Path
from typing import Any

from lib.assembly import Assembly
from ..CHM13v2.assembly import CHM13v2
from ..HSV1.assembly import HSV1
from ..IAV.assembly import IAV


class zdhs(Assembly):
    def __init__(
            self, fasta: Path, gff3: Path, hsv1: HSV1, iav: IAV, chm13v2: CHM13v2
    ):
        super().__init__()
        self.fasta = fasta
        self.gff3 = gff3
        self.hsv1 = hsv1
        self.iav = iav
        self.chm13v2 = chm13v2

    @property
    def name(self) -> str:
        return "zdhs"

    @property
    def organisms(self) -> frozenset[str]:
        return frozenset([*self.hsv1.organisms, *self.iav.organisms, *self.chm13v2.organisms])

    def seqsizes(self) -> dict[str, int]:
        return self.hsv1.seqsizes() | self.iav.seqsizes() | self.chm13v2.seqsizes()

    def seqlabels(self) -> dict[str, str]:
        return self.hsv1.seqlabels() | self.iav.seqlabels() | self.chm13v2.seqlabels()

    def annotations(self) -> frozenset[str]:
        return frozenset([
            "REDIportal", "RepeatMasker", "RepeatMasker classification", "GFF3", "FASTA",
        ])

    def files(self, annotation: str) -> tuple[Path, ...]:
        match annotation:
            case "REDIportal" | "RepeatMasker" | "RepeatMasker classification":
                return self.chm13v2.files(annotation)
            case "GFF3":
                return (self.gff3,)
            case "FASTA":
                return (self.fasta,)
            case _:
                raise ValueError(f"Unknown annotation: {annotation}")

    def load(self, annotation: str) -> Any:
        match annotation:
            case "REDIportal" | "RepeatMasker" | "RepeatMasker classification":
                return self.chm13v2.load(annotation)
            case "GFF3" | "FASTA":
                raise NotImplementedError(f"Loading {annotation} is unsupported")
            case _:
                raise ValueError(f"Unknown annotation: {annotation}")
