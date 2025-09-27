import bisect
from collections.abc import Hashable
from typing import Optional, Iterable, Literal

from biobit.core.loc import Interval


class BoundariesIndex[T: Hashable]:
    def __init__(self, boundaries: Iterable[tuple[int, T]]):
        collapsed = set()
        for b in boundaries:
            collapsed |= set(b)

        self.boundaries = []
        self.data = []
        for pos, data in sorted(boundaries):
            self.boundaries.append(pos)
            self.data.append(data)

    def closest(self, pos: int, side: Literal['left', 'right']) -> tuple[int, T]:
        assert self.boundaries[0] <= pos <= self.boundaries[-1], (pos, self.boundaries[0], self.boundaries[-1])
        ind = bisect.bisect_left(self.boundaries, pos)  # self.boundaries[ind] >= pos
        match side:
            case 'left':
                while self.boundaries[ind - 1] == pos:
                    ind -= 1
                neighbor = self.boundaries[ind - 1]
                return neighbor, self.data[ind - 1]
            case 'right':
                while self.boundaries[ind] == pos:
                    ind += 1
                neighbor = self.boundaries[ind]
                return neighbor, self.data[ind]
            case _:
                raise ValueError(side)

    def window(self, pos: int, side: Literal['left', 'right'], maxsize: Optional[int] = None) -> tuple[Interval, T]:
        neighbor, data = self.closest(pos, side)
        match side:
            case 'left':
                if maxsize:
                    neighbor = max(neighbor, pos - maxsize)
                window = Interval(neighbor, pos)
            case 'right':
                if maxsize:
                    neighbor = min(neighbor, pos + maxsize)
                window = Interval(pos, neighbor)
            case _:
                raise ValueError(side)
        return window, data

    def __contains__(self, item):
        ind = bisect.bisect_left(self.boundaries, item)
        return self.boundaries[ind] == item

    # def remove(self, positions: Iterable[int]):
    #     self._boundaries_set -= set(positions)
    #     self.boundaries = sorted(self._boundaries_set)


class TranscriptionBoundaries(BoundariesIndex[Literal['tss', 'tes', 'donor', 'acceptor', 'seq-start', 'seq-end']]):
    ...
