from pathlib import Path
from typing import Any

from biobit.deprecated.repmasker import RepmaskerClassification
from pybedtools import BedTool

from lib.assembly import Assembly, HasRepeatMasker, GencodeLiftoffAnnotome, HasRefSeqAnnotome, \
    HasGencodeLiftoffAnnotome, RefSeqAnnotome
from lib.pkl import PklData
from . import seqid


class CHM13v2(Assembly, HasGencodeLiftoffAnnotome, HasRefSeqAnnotome, HasRepeatMasker):
    def __init__(
            self,
            fasta: Path,
            gencode_liftoff_gff3: Path, gencode_liftoff_annotome: PklData[GencodeLiftoffAnnotome],
            refseq_gff3: Path, refseq_annotome: PklData[RefSeqAnnotome],
            rediportal: Path,
            repmasker: Path, repmasker_classification: Path
    ):
        super().__init__()
        self._fasta = fasta
        self._gencode_liftoff_gff3 = gencode_liftoff_gff3
        self._gencode_liftoff_annotome = gencode_liftoff_annotome
        self._refseq_gff3 = refseq_gff3
        self._refseq_annotome = refseq_annotome
        self._rediportal = rediportal
        self._repmasker = repmasker
        self._repmasker_classification = repmasker_classification

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
                return (self._rediportal,)
            case "RepeatMasker":
                return (self._repmasker,)
            case "RepeatMasker classification":
                return (self._repmasker_classification,)
            case "GENCODE-Liftoff":
                return (self._gencode_liftoff_gff3,)
            case "RefSeq":
                return (self._refseq_gff3,)
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
                return RepmaskerClassification(self._repmasker_classification)
            case "GENCODE-Liftoff":
                return self._gencode_liftoff_annotome.load()
            case "RefSeq":
                return self._refseq_annotome.load()
            case "FASTA":
                raise NotImplementedError(f"Loading {annotation} is not implemented yet.")
            case _:
                raise ValueError(f"Unknown annotation: {annotation}")

    def gencode_liftoff(self) -> GencodeLiftoffAnnotome:
        return self.load("GENCODE-Liftoff")

    def refseq(self) -> RefSeqAnnotome:
        return self.load("RefSeq")

    def repeat_masker(self) -> BedTool:
        return self.load("RepeatMasker")

    def repeat_masker_classification(self) -> RepmaskerClassification:
        return self.load("RepeatMasker classification")
