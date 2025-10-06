from collections import defaultdict
from dataclasses import dataclass
from typing import Iterator

from biobit.core.loc import IntoLocus, Orientation, Interval, Locus
from biobit.toolkit.repeto.repeats import InvRepeat


@dataclass(frozen=True, slots=True)
class RepetoGroup:
    contig: str
    orientation: Orientation
    bsegment: Interval

    rois: list[Interval]

    def __post_init__(self):
        expstart, expend = min(x.start for x in self.rois), max(x.end for x in self.rois)
        assert (expstart, expend) == self.bsegment, ((expstart, expend), self.bsegment, self.rois)

    def roi_pairs(self, maxdist: int) -> Iterator[tuple[Interval, Interval]]:
        for i in range(len(self.rois)):
            for j in range(i + 1, len(self.rois)):
                if self.rois[j].start - self.rois[i].end <= maxdist:
                    yield self.rois[i], self.rois[j]


def group(
        rois: list[IntoLocus], insulators: list[IntoLocus], connectors: list[IntoLocus], maxdist: int
) -> list[RepetoGroup]:
    # 1. Group based on contig and orientation
    _rois, _insulators, _connectors = defaultdict(list), defaultdict(list), defaultdict(list)
    for regions, saveto in (rois, _rois), (insulators, _insulators), (connectors, _connectors):
        for region in regions:
            if isinstance(region, tuple):
                region = Locus(*region)
            saveto[region.contig, region.orientation].append(region.interval)

    # 2. Group based on distance
    skipped_insulators = 0
    total_insulators = len(insulators)

    result = []
    for (seqid, orientation), rois in _rois.items():
        insulators = _insulators.get((seqid, orientation), [])
        connectors = _connectors.get((seqid, orientation), [])

        rois = Interval.merge(rois + connectors)
        insulators = Interval.merge(insulators)

        # Sanity check
        assert all(nxt.start > prv.end for nxt, prv in zip(rois[1:], rois[:-1]))
        assert all(nxt.start > prv.end for nxt, prv in zip(insulators[1:], insulators[:-1]))

        cache = [rois[0]]
        start, end = rois[0].start, rois[0].end

        total_insulators += len(insulators)
        insulate = insulators[0] if insulators else None
        insulators = insulators[1:]

        for roi in rois[1:]:
            distance = roi.start - end

            if distance > maxdist:
                # ROIs are far away - no grouping needed
                result.append(RepetoGroup(seqid, orientation, Interval(start, end), cache))
                cache = [roi]
                start, end = roi.start, roi.end
            else:
                # ROIs are close enough, maybe we need to group them
                if insulate:
                    if insulate.intersects(roi) or insulate.intersects((start, end)):
                        skipped_insulators += 1
                        end = roi.end
                        cache.append(roi)
                    if end <= insulate.start and roi.start >= insulate.end:
                        # Insulator in between - reset the cache
                        result.append(RepetoGroup(seqid, orientation, Interval(start, end), cache))
                        cache = [roi]
                        start, end = roi.start, roi.end
                    else:
                        # No insulator in between - just extend the cache
                        end = roi.end
                        cache.append(roi)
                else:
                    # No insulator in between - just extend the cache
                    end = roi.end
                    cache.append(roi)

            # Fast-forward insulators if needed
            while insulators and insulate and roi.start >= insulate.end:
                insulate = insulators[0]
                insulators = insulators[1:]

        result.append(RepetoGroup(seqid, orientation, Interval(start, end), cache))

    # Drop connectors from the result
    excluding_connectors = []
    for group in result:
        # Any overlap with connectors?
        if (group.contig, group.orientation) not in _connectors:
            excluding_connectors.append(group)
            continue

        connectors = Interval.overlap([group.bsegment], _connectors[group.contig, group.orientation])
        if not connectors:
            excluding_connectors.append(group)
            continue

        # Brute force recalculate the ROIs
        rois = Interval.overlap(_rois[group.contig, group.orientation], [group.bsegment])
        rois = Interval.merge(rois)
        if rois:
            object.__setattr__(group, 'rois', rois)
            excluding_connectors.append(group)

    result = excluding_connectors

    if skipped_insulators > 0:
        print(f"Skipped insulators: {skipped_insulators} ({skipped_insulators / total_insulators:.1%})")

    return result


@dataclass(frozen=True, slots=True)
class Partition:
    ind: str
    contig: str
    orientation: Orientation

    invrep: list[InvRepeat]
    peaks: list[Interval]

    def envelope(self) -> Interval:
        allsegments = [self.peaks]
        allsegments.extend(x.seqranges() for x in self.invrep)

        start = min(x.start for segments in allsegments for x in segments)
        end = max(x.end for segments in allsegments for x in segments)
        return Interval(start, end)
