from dataclasses import dataclass
from typing import Literal

from biobit.core.loc import Strand, Interval


@dataclass(frozen=True, slots=True)
class NormBin:
    # Unique identifier of the bin and its type
    ind: str
    type: Literal['exonic', 'intronic', 'intergenic', 'vRNA', 'MT']

    # Location of the partition
    seqid: str
    strand: Strand
    intervals: tuple[Interval, ...]
