from attrs import define
from biobit.core.loc import Strand, Interval


@define(frozen=True)
class IntronCore:
    # Donor and acceptor exons
    donor: Interval
    acceptor: Interval
    # Inside long introns there might multiple intronic "pieces" gapped by small RNAs, nested genes, etc.
    intron: tuple[Interval, ...]

    def __attrs_post_init__(self):
        allsegments = sorted([self.donor, *self.intron, self.acceptor], key=lambda x: x.start)
        for prv, nxt in zip(allsegments[:-1], allsegments[1:]):
            assert prv.start < prv.end <= nxt.start < nxt.end, f"{prv} overlaps {nxt}"


@define(frozen=True)
class RNACore:
    gid: str
    seqid: str
    strand: Strand

    # Ensembl canonical transcript's exons trimmed at the 5` and 3` ends to match 'shortest' TSS and TES
    exons: tuple[Interval, ...]
    # Canonical transcript's introns gapped by small RNAs, nested genes, etc.
    introns: tuple[IntronCore, ...] | None

    # Read-through region (downstream of all annotated TESs)
    read_through: Interval
    # Read-in region (upstream of all annotated TSSs)
    read_in: Interval
    # Divergent transcription
    divergent: Interval

    def truncate(self, max_read_through: int, max_read_in: int, max_divergent: int) -> 'RNACore':
        match self.strand:
            case "+":
                read_through = Interval(
                    self.read_through.start, min(self.read_through.end, self.read_through.start + max_read_through)
                )
                readin = Interval(max(self.read_in.start, self.read_in.end - max_read_in), self.read_in.end)
                divergent = Interval(max(self.divergent.start, self.divergent.end - max_divergent), self.divergent.end)
            case "-":
                read_through = Interval(
                    max(self.read_through.start, self.read_through.end - max_read_through), self.read_through.end
                )
                readin = Interval(self.read_in.start, min(self.read_in.end, self.read_in.start + max_read_in))
                divergent = Interval(self.divergent.start,
                                     min(self.divergent.end, self.divergent.start + max_divergent))
            case _:
                raise ValueError(f"Invalid strand {self.strand}")
        return RNACore(
            self.gid, self.seqid, self.strand, self.exons, self.introns, read_through, readin, divergent
        )
