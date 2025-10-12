from abc import ABCMeta, abstractmethod
from pathlib import Path
from typing import Any, runtime_checkable, Protocol

from biobit.deprecated.repmasker import RepmaskerClassification
from pybedtools import BedTool

from .annotation.gencode import GencodeAnnotome, GencodeLiftoffAnnotome
from .annotation.refseq import RefSeqAnnotome


class Assembly(metaclass=ABCMeta):
    """A formal contract for a genomic assembly."""

    @property
    @abstractmethod
    def name(self) -> str:
        """A unique name for the assembly (e.g., 'CHM13v2')."""
        ...

    @property
    @abstractmethod
    def organisms(self) -> frozenset[str]:
        """The set of organisms covered by the assembly (e.g., {'Homo sapiens', 'HSV-1'})."""
        ...

    @abstractmethod
    def seqsizes(self) -> dict[str, int]:
        """Returns a mapping of sequence IDs to their lengths."""
        ...

    def seqlabels(self) -> dict[str, str]:
        """Returns a mapping of sequence names to their human-readable labels whenever possible."""
        return {k: k for k in self.seqsizes()}

    @abstractmethod
    def annotations(self) -> frozenset[str]:
        """Returns the set of available annotation keys."""
        ...

    @abstractmethod
    def files(self, annotation: str) -> tuple[Path, ...]:
        """
        Gets the raw file path(s) for a given annotation key.

        Raises KeyError if the key is not found.
        """
        ...

    @abstractmethod
    def load(self, annotation: str) -> Any:
        """
        Loads and parses an annotation, returning a data structure.

        Raises KeyError if the key is not found.
        """
        ...


@runtime_checkable
class HasRepeatMasker(Protocol):
    def repeat_masker(self) -> BedTool: ...

    def repeat_masker_classes(self) -> RepmaskerClassification: ...


@runtime_checkable
class HasRefSeqAnnotome(Protocol):
    def refseq(self) -> RefSeqAnnotome: ...


@runtime_checkable
class HasGencodeAnnotome(Protocol):
    def gencode(self) -> GencodeAnnotome: ...


@runtime_checkable
class HasGencodeLiftoffAnnotome(Protocol):
    def gencode_liftoff(self) -> GencodeLiftoffAnnotome: ...
