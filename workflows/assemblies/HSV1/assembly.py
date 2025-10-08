from pathlib import Path
from typing import Any, Literal

from biobit.core.loc import Interval

from .resources import expression
from ..resources import Assembly


class HSV1(Assembly):
    def __init__(self, fasta: Path, annotation_gff3: Path):
        super().__init__()
        self.fasta = fasta
        self.annotation_gff3 = annotation_gff3

    @property
    def name(self) -> str:
        return "HSV1"

    @property
    def organisms(self) -> frozenset[str]:
        return frozenset(["Herpes simplex virus 1"])

    def seqsizes(self) -> dict[str, int]:
        return {"NC_001806.2": 152222}

    def annotations(self) -> frozenset[str]:
        return frozenset(["GFF3", "FASTA"])

    def files(self, annotation: str) -> tuple[Path, ...]:
        match annotation:
            case "GFF3":
                return (self.annotation_gff3,)
            case "FASTA":
                return (self.fasta,)
            case _:
                raise ValueError(f"Unknown annotation: {annotation}")

    def load(self, annotation: str) -> Any:
        raise NotImplementedError("Loading of annotations is not implemented.")

    def expression_stage(self, rna: str) -> Literal["Latency", "Immediate Early", "Early", "Late", "Uncharacterized"]:
        return expression.STAGES.get(rna, "Uncharacterized")

    def nagnag_spliced_rnas(self) -> set[str]:
        return {"TRL2-iso2", "UL36.6-iso2", "IRL2-iso2"}

    def terminal_repeats(self) -> tuple[Interval, Interval]:
        return Interval(0, 9213), Interval(145589, 152222)
