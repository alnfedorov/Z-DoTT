from pathlib import Path
from typing import Any

from pybedtools import BedTool

from lib.assembly import Assembly, HasGencodeAnnotome, GencodeAnnotome, HasRefSeqAnnotome, RefSeqAnnotome
from lib.pkl import PklData
from . import seqid


class GRCm39(Assembly, HasGencodeAnnotome, HasRefSeqAnnotome):
    def __init__(
            self,
            fasta: Path,
            gencode_gff3: Path, gencode_annotome: PklData[GencodeAnnotome],
            refseq_gff3: Path, refseq_annotome: PklData[RefSeqAnnotome],
            rediportal: Path,
            repmasker: Path, repmasker_classification: Path
    ):
        super().__init__()
        self._fasta = fasta
        self._gencode = gencode_gff3
        self._gencode_annotome = gencode_annotome
        self._refseq = refseq_gff3
        self._refseq_annotome = refseq_annotome
        self._rediportal = rediportal
        self._repmasker = repmasker
        self._repmasker_classification = repmasker_classification

    @property
    def name(self) -> str:
        return "GRCm39"

    @property
    def organisms(self) -> frozenset[str]:
        return frozenset(["Mus musculus"])

    def seqsizes(self) -> dict[str, int]:
        return seqid.sizes()

    def annotations(self) -> frozenset[str]:
        return frozenset([
            "REDIportal", "RepeatMasker", "RepeatMasker classification", "GENCODE", "RefSeq", "FASTA",
        ])

    def files(self, annotation: str) -> tuple[Path, ...]:
        match annotation:
            case "REDIportal":
                return (self._rediportal,)
            case "RepeatMasker":
                return (self._repmasker,)
            case "RepeatMasker classification":
                return (self._repmasker_classification,)
            case "GENCODE":
                return (self._gencode,)
            case "RefSeq":
                return (self._refseq,)
            case "FASTA":
                return (self._fasta,)
            case _:
                raise ValueError(f"Unknown annotation: {annotation}")

    def load(self, annotation: str) -> Any:
        match annotation:
            case "REDIportal":
                return BedTool(self._rediportal.as_posix())
            case "RepeatMasker":
                return BedTool(self._repmasker.as_posix())
            case "RepeatMasker classification":
                return self._repmasker_classification
            case "GENCODE":
                return self._gencode_annotome.load()
            case "RefSeq":
                return self._refseq_annotome.load()
            case "FASTA":
                raise NotImplementedError(f"Loading {annotation} is unsupported")
            case _:
                raise ValueError(f"Unknown annotation: {annotation}")

    def gencode(self) -> GencodeAnnotome:
        return self.load("GENCODE")

    def refseq(self) -> RefSeqAnnotome:
        return self.load("RefSeq")
