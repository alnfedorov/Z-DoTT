from pathlib import Path
from typing import Any

from ..resources import Assembly


class IAV(Assembly):
    def __init__(self, fasta: Path, annotation_gff3: Path):
        super().__init__()
        self._fasta = fasta
        self._gff3 = annotation_gff3

    @property
    def name(self) -> str:
        return "IAV"

    @property
    def organisms(self) -> frozenset[str]:
        return frozenset(["Influenza A virus"])

    def seqsizes(self) -> dict[str, int]:
        return {
            "NC_002023.1": 2341,
            "NC_002022.1": 2233,
            "NC_002021.1": 2341,
            "NC_002020.1": 890,
            "NC_002019.1": 1565,
            "NC_002018.1": 1413,
            "NC_002017.1": 1778,
            "NC_002016.1": 1027,
        }

    def annotations(self) -> frozenset[str]:
        return frozenset(["GFF3", "FASTA"])

    def files(self, annotation: str) -> tuple[Path, ...]:
        match annotation:
            case "GFF3":
                return (self._gff3,)
            case "FASTA":
                return (self._fasta,)
            case _:
                raise ValueError(f"Unknown annotation: {annotation}")

    def load(self, annotation: str) -> Any:
        raise NotImplementedError("Loading of annotations is not implemented.")
