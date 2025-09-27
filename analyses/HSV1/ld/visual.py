from collections import defaultdict
from typing import Callable, List, Optional, Tuple

import numpy as np
import numpy.typing as npt
from pybedtools import Interval


def simplify(values: npt.NDArray[np.float32],
             merge: Callable[[List[float]], float] = np.mean,
             abseps: Optional[float] = None,
             releps: Optional[float] = None,
             ndigits: int = 3) -> List[Tuple[int, int, float]]:
    assert abseps or releps

    result = []
    start, cache = 0, [values[0]]
    for coord, nextval in enumerate(values[1:]):
        # account for the skipped first element
        coord += 1

        curval = merge(cache)
        diff = abs(nextval - curval)
        if (abseps and diff > abseps) or (releps and curval > 0 and abs(diff / curval) > releps):
            result.append((start, coord, round(curval, ndigits)))
            start, cache = coord, [nextval]
        else:
            cache.append(nextval)

    result.append((start, len(values), round(merge(cache), ndigits)))
    return result


def assort_genes_by_height(intervals: List[Interval],
                           limits: Tuple[float, float] = (0.8, 0.2),
                           offset: int = 300,
                           step: Optional[float] = None) -> List[Interval]:
    if len(intervals) == 0:
        return []

    intervals = sorted(intervals, key=lambda x: x.start)
    # group by contig
    contigs = defaultdict(list)
    for i in intervals:
        contigs[i.chrom].append(i)

    levels = {}
    for contig, intervals in contigs.items():
        intervals = sorted(intervals, key=lambda x: (x.start, -x.length))
        contigs[contig] = intervals

        # Greedy fill intervals level-by-level
        curlvl = 0
        left = list(range(len(intervals)))
        while True:
            if len(left) == 0:
                break

            curind, nextiter = left[0], []
            for ind in left[1:]:
                assert intervals[ind].start >= intervals[curind].start
                # Can we pick the next interval?
                if intervals[ind].start - intervals[curind].end > offset:
                    levels[(curind, contig)] = curlvl
                    curind = ind
                else:
                    nextiter.append(ind)
            levels[(curind, contig)] = curlvl

            curlvl += 1
            left = nextiter

    maxlvl = max(levels.values())
    if step is None:
        if maxlvl == 0:
            step = 0
        else:
            step = abs(limits[1] - limits[0]) / maxlvl

    result = []
    for intervals in contigs.values():
        for ind, x in enumerate(intervals):
            lvl = levels[ind, x.chrom] * step
            if limits[0] > limits[1]:
                lvl = limits[0] - lvl
            else:
                lvl = limits[0] + lvl
            # assert round(min(limits), 3) <= round(lvl, 3) <= round(max(limits), 3)
            x.score = str(lvl)
            result.append(x)
    return result
