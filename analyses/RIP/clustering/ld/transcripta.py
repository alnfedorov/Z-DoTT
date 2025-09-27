from collections import defaultdict
from dataclasses import dataclass, field

from biobit.core.loc import Interval, Strand
from biobit.core.loc.mapping import ChainMap
from biobit.toolkit.annotome import Annotome
from intervaltree import IntervalTree

from stories import annotation


@dataclass(frozen=True, slots=True)
class Transcript:
    ind: str
    contig: str
    strand: Strand

    exons: list[Interval]
    mapping: ChainMap = field(init=False)

    def __post_init__(self):
        object.__setattr__(self, "exons", sorted(self.exons, key=lambda x: x.start))
        object.__setattr__(self, "mapping", ChainMap(Interval.merge(self.exons)))

    def bbox(self) -> Interval:
        return Interval(self.exons[0].start, self.exons[-1].end)

    def map(self, segment: Interval) -> Interval | None:
        return self.mapping.map_interval(segment)

    def __hash__(self):
        return hash((self.ind, self.contig, self.strand))


def parse(assembly) -> dict[tuple[str, Strand], IntervalTree]:
    index = defaultdict(IntervalTree)

    gencode: Annotome = assembly.gencode.load()
    for rna in gencode.rnas.values():
        if annotation.filters.is_primary(rna) and rna.attrs.type == "protein_coding":
            transcript = Transcript(rna.ind, rna.loc.seqid, rna.loc.strand, rna.exons)
            bbox = transcript.bbox()
            index[rna.loc.seqid, rna.loc.strand].addi(bbox.start, bbox.end, data=transcript)

    return index
