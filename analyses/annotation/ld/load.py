import pickle
from functools import lru_cache

from biobit.core.loc import Strand, Orientation, Interval
from biobit.toolkit.annotome import Annotome

from . import paths
from .boundaries import TranscriptionBoundaries
from .rna_core import RNACore


@lru_cache(maxsize=None)
def gencode[AttrGene, AttrRNA, AttrCDS](assembly: str) -> Annotome[AttrGene, AttrRNA, AttrCDS]:
    with open(paths.gencode.pkl[assembly], 'rb') as stream:
        return pickle.load(stream)


@lru_cache(maxsize=None)
def rna_cores(assembly: str) -> dict[tuple[str, Strand], list[RNACore]]:
    with open(paths.rna_cores.pkl[assembly], 'rb') as stream:
        return pickle.load(stream)


@lru_cache(maxsize=None)
def transcription_boundaries(assembly: str) -> dict[tuple[str, Strand], TranscriptionBoundaries]:
    with open(paths.transcription_boundaries.pkl[assembly], 'rb') as stream:
        return pickle.load(stream)


@lru_cache(maxsize=None)
def resolved_annotation(assembly: str) -> dict[str, dict[tuple[str, Orientation], list[Interval]]]:
    with open(paths.resolved_annotation.pkl[assembly], 'rb') as stream:
        return pickle.load(stream)
