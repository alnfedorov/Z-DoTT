from bisect import bisect_left
from collections import defaultdict
from itertools import chain
from typing import Optional

import numpy as np
import numpy.typing as npt
import pyBigWig
from attr import dataclass, field
from biobit.core.loc import Strand, Interval
from biobit.toolkit import nfcore
from biobit.toolkit.repeto.repeats import InvRepeat, InvSegment
from intervaltree import IntervalTree

from stories.RIP import pcalling


class ExperimentTracks:
    def __init__(self, cmp: pcalling.Config):
        signal, control = {Strand.Forward: [], Strand.Reverse: []}, {Strand.Forward: [], Strand.Reverse: []}
        for exps, saveto in (cmp.signal, signal), (cmp.control, control):
            for exp in exps:
                fwd, rev = nfcore.rnaseq.extract.bigwig(exp)

                assert fwd.is_file(), fwd
                saveto[Strand.Forward].append(fwd)

                assert rev.is_file(), rev
                saveto[Strand.Reverse].append(rev)

        self._contig: Optional[str] = None
        self._ctglen: Optional[int] = None
        self._strand: Optional[Strand] = None

        self._signal_path = signal
        self._signal_values: Optional[npt.NDArray[np.float32]] = None

        self._control_path = control
        self._control_values: Optional[npt.NDArray[np.float32]] = None

    def open(self, contig: str, ctglen: int, strand: Strand) -> 'ExperimentTracks':
        self._contig = contig
        self._ctglen = ctglen
        self._strand = strand

        # Load the whole contig into memory and calculate the average signal/control across all experiments
        results = []
        for paths in self._control_path[strand], self._signal_path[strand]:
            values = np.zeros(ctglen, dtype=np.float32)
            for p in paths:
                with pyBigWig.open(p.as_posix()) as bw:
                    fetched = bw.values(contig, 0, ctglen, numpy=True)
                    fetched[np.isnan(fetched)] = 0
                    values += fetched
                del fetched
            values /= len(paths)
            results.append(values)
        self._control_values, self._signal_values = results
        return self

    def signal(self, start: int, end: int):
        return self._signal_values[start:end]

    def control(self, start: int, end: int):
        return self._control_values[start:end]

    def score(self, start: int, end: int, insulators: list[tuple[int, int]], peaks: IntervalTree) -> 'ScoreState':
        diff = self.signal(start, end) - self.control(start, end)
        diff -= diff.min()

        insulators = sorted(x for y in insulators for x in y)
        return ScoreState(start, end, diff, insulators, peaks)

    def __getstate__(self):
        return self._contig, self._ctglen, self._strand, self._signal_path, self._control_path

    def __setstate__(self, state):
        self._contig, self._ctglen, self._strand, self._signal_path, self._control_path = state
        if self._contig and self._strand:
            self.open(self._contig, self._ctglen, self._strand)


@dataclass(frozen=True)
class ScoreState:
    start: int
    end: int
    scores: npt.NDArray[np.float32]
    insulators: list[int]
    peaks: IntervalTree
    max_overlap: int = field(init=False)

    def __attrs_post_init__(self):
        total_length = sum(x.length() for x in self.peaks)
        object.__setattr__(self, "max_overlap", total_length)

    def _score(self, left: Interval, right: Interval):
        lscore = self.scores[left.start - self.start: left.end - self.start]
        rscore = self.scores[right.start - self.start: right.end - self.start][::-1]
        return np.minimum(lscore, rscore)

    def score(self, ir: InvRepeat) -> float:
        # Raw score is the total signal across dsRNA arms
        score = sum(self._score(segment.left, segment.right).sum() for segment in ir.segments)
        score = float(score)

        weight = 1.0
        left, right = ir.left_brange(), ir.right_brange()

        # Downscale the score based on the overlap between dsRNA arms and sample peaks
        overlap = 0
        for segment in left, right:
            for ov in self.peaks.overlap(segment.start, segment.end):
                overlap += ov.overlap_size(segment.start, segment.end)
        weight *= overlap / self.max_overlap

        # Downscale the score based on the distance between the two arms
        distance = right.start - left.end
        if distance > 1_000:
            weight /= (distance / 1_000)

        # Downscale the score if there is an insulator between the two arms
        lind = bisect_left(self.insulators, left.end)
        rind = bisect_left(self.insulators, right.start)
        if lind != rind:
            weight /= 2 ** min(20, abs(lind - rind))

        score *= weight

        return score

    def resolve(self, ir: InvRepeat) -> float:
        score = 0
        for segment in ir.segments:
            left, right = segment.left, segment.right
            explained = self._score(left, right)

            score += float(explained.sum())

            for segment, sub in (left, explained), (right, explained[::-1]):
                self.scores[segment.start - self.start: segment.end - self.start] = np.maximum(
                    0, self.scores[segment.start - self.start: segment.end - self.start] - sub
                )
        return score


def from_dsRNA_coordinates_to_global(rna: InvRepeat, intervals: list[Interval]) -> list[InvRepeat]:
    result = []

    start, end, ind = 0, rna.segments[0].left.len(), 0
    for iv in sorted(intervals):
        # Fast-forward to the segment that contains the overlap
        while end <= iv.start:
            start, end, ind = end, end + rna.segments[ind + 1].left.len(), ind + 1
        assert start <= iv.start < end

        cache = []

        # Crop the first segment
        loffset, roffset = iv.start - start, end - min(iv.end, end)
        left = Interval(rna.segments[ind].left.start + loffset, rna.segments[ind].left.end - roffset)
        right = Interval(rna.segments[ind].right.start + roffset, rna.segments[ind].right.end - loffset)
        cache.append(InvSegment(left, right))

        # Stop if the overlap is on the last segment or if the next segment is not overlapping
        if ind + 1 == len(rna.segments) or iv.end < end:
            assert sum(x.left.len() for x in cache) == iv.len(), (iv, cache)
            ir = InvRepeat(cache)
            result.append(ir)
            continue

        # Otherwise, add all segments that are fully included in the overlap
        start, end, ind = end, end + rna.segments[ind + 1].left.len(), ind + 1
        while end <= iv.end and (ind + 1) < len(rna.segments):
            assert iv.start <= start < end <= iv.end
            cache.append(rna.segments[ind])
            start, end, ind = end, end + rna.segments[ind + 1].left.len(), ind + 1

        # Crop the last segment if needed
        if start < iv.end <= end:
            offset = end - iv.end
            left = Interval(rna.segments[ind].left.start, rna.segments[ind].left.end - offset)
            right = Interval(rna.segments[ind].right.start + offset, rna.segments[ind].right.end)
            cache.append(InvSegment(left, right))

        # Save the results
        assert sum(x.left.len() for x in cache) == iv.len(), (iv, cache)
        ir = InvRepeat(cache)
        result.append(ir)
    return result


def filter_segments(rna: InvRepeat, solutions: list[InvRepeat], min_samples: int):
    solutions = sorted(segment.left for solution in solutions for segment in solution.segments)

    # Map solutions to dsRNA coordinates
    # Suffice to map only the left arm, the right arm is always the same
    mapped, pos = [], 0
    for segment in rna.segments:
        for solution in solutions:
            if (ov := segment.left.intersection(solution)) is not None:
                mapped.append(Interval(ov.start - segment.left.start + pos, ov.end - segment.left.start + pos))
        pos += segment.left.len()
    mapped = sorted(mapped)

    # Select regions that are support by at least X samples
    # Very slow algorithm, but it's not a bottleneck
    replicated = defaultdict(list)
    boundaries = sorted(set(chain.from_iterable((x.start, x.end) for x in mapped)))
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        interval = Interval(start, end)
        count = sum(x.intersects(interval) for x in mapped)
        if count >= min_samples:
            replicated[count].append(interval)

    if not replicated:
        return [], []

    # # Merge-Map each region back to dsRNA coordinates
    # backmapped, tags = [], []
    # for count, intervals in replicated.items():
    #     bckmp = from_dsRNA_coordinates_to_global(rna, Interval.merge(intervals))
    #     bckmp = sorted([segment for ir in bckmp for segment in ir.segments], key=lambda x: x.left.start)
    #
    #     backmapped.append(InvRepeat(bckmp))
    #     tags.append(count)

    # Merge the intervals and map back to dsRNA coordinates
    intervals = [it for intervals in replicated.values() for it in intervals]
    backmapped = from_dsRNA_coordinates_to_global(rna, Interval.merge(intervals))

    segments = sorted([segment for ir in backmapped for segment in ir.segments], key=lambda x: x.left.start)
    backmapped = [InvRepeat(segments)]
    tags = [f"N>={min(replicated.keys())}"]

    return backmapped, tags
