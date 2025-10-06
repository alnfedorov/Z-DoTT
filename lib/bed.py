from collections import defaultdict
from pathlib import Path
from typing import Iterable, Callable, Any, Optional

import pybedtools
from biobit.core.loc import Interval, Orientation, IntoOrientation


def tbindex(bed: pybedtools.BedTool, saveto: Path):
    assert saveto.suffixes[-2:] == [".bed", ".gz"]

    saveto = saveto.with_suffix("")
    bed.saveas(saveto).tabix(in_place=True, force=True, is_sorted=True)
    saveto.unlink()


def group(
        intervals: Iterable[Interval],
        key: Callable[[Interval], Any] = lambda x: (x.chrom, x.strand),
        value: Callable[[Interval], Any] = lambda x: (x.start, x.end)
) -> dict[Any, list[Any]]:
    group = defaultdict(list)
    for it in intervals:
        group[key(it)].append(value(it))
    return dict(group)


class blocks:
    @staticmethod
    def make(
            seqid: str, intervals: Iterable[Interval], orientation: IntoOrientation, name: str, color: str = "0,0,0"
    ) -> pybedtools.Interval:
        intervals = sorted(intervals, key=lambda x: x.start)

        assert len(intervals) >= 1
        for ind, (prv, nxt) in enumerate(zip(intervals[:-1], intervals[1:])):
            assert prv.start < prv.end <= nxt.start < nxt.end, f"{prv} overlaps {nxt}"

        blockCount = str(len(intervals))
        blockSizes = ",".join(map(str, [s.end - s.start for s in intervals]))
        blockStarts = ",".join(map(str, [s.start - intervals[0].start for s in intervals]))

        start, end = intervals[0].start, intervals[-1].end
        thickStart, thickEnd = str(start), str(end)
        return pybedtools.Interval(
            seqid, start, end, name, strand=Orientation(orientation).symbol(),
            otherfields=[thickStart, thickEnd, color, blockCount, blockSizes, blockStarts]
        )

    @staticmethod
    def from_iterable(
            iterator: Iterable[tuple[str, Iterable[Interval], IntoOrientation, str, Optional[str]]],
            sort: bool = True
    ) -> pybedtools.BedTool:
        bed = []
        for seqid, intervals, strand, name, color in iterator:
            if color is None:
                color = "0,0,0"
            bed.append(blocks.make(seqid, intervals, strand, name, color))
        bed = pybedtools.BedTool(bed)
        if sort:
            bed = bed.sort()
        return bed


def parse_coords(coords: str, **kwargs) -> pybedtools.Interval:
    coords = coords.replace(',', '')
    if ":" in coords:
        seqid, extra = coords.split(':')
        start, end = extra.split('-')
    else:
        seqid, start, end = coords.split()
    return pybedtools.Interval(seqid, int(start), int(end), **kwargs)


def merge_stranded(bed: list[pybedtools.BedTool], **kwargs) -> pybedtools.BedTool:
    assert len(bed) >= 1
    if len(bed) == 1:
        bed = bed[0]
    else:
        bed = bed[0].cat(*bed[1:], postmerge=False)

    if len(bed) == 0:
        return bed

    bed = bed.sort().merge(s=True, c=6, o="distinct", **kwargs)
    bed = pybedtools.BedTool([pybedtools.Interval(x.chrom, x.start, x.end, strand=x.fields[3]) for x in bed])
    return bed.sort()
