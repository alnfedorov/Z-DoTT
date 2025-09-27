import pickle
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
from biobit.toolkit import seqproj, reaper
from biobit.toolkit.libnorm import MedianOfRatiosNormalization

PKL = Path(__file__).with_suffix(".pkl")


@dataclass(frozen=True)
class ReaperFiles:
    signal: Path
    control: Path
    modeled: Path
    raw_peaks: Path
    filtered_peaks: Path

    def __iter__(self):
        return iter((self.signal, self.control, self.modeled, self.raw_peaks, self.filtered_peaks))


@dataclass(frozen=True)
class Config:
    ind: str

    project: str
    host: str
    organism: frozenset[str]
    assembly: str

    signal: tuple[seqproj.Experiment, ...]
    control: tuple[seqproj.Experiment, ...]

    model: reaper.model.RNAPileup
    pcalling: reaper.pcalling.ByCutoff
    nms: reaper.postfilter.NMS

    root: Path
    reaper: ReaperFiles = field(init=False)

    def __post_init__(self):
        title = self.ind.replace("/", "-").replace(" ", "-")

        signal = self.root / Path("signal") / f"{title}.signal.bed.gz"
        control = self.root / Path("control") / f"{title}.control.bed.gz"
        modeled = self.root / Path("model") / f"{title}.modeled.bed.gz"
        raw_peaks = self.root / Path("raw-peaks") / f"{title}.raw.bed.gz"
        filtered_peaks = self.root / Path("filtered-peaks") / f"{title}.filtered.bed.gz"

        object.__setattr__(self, "reaper", ReaperFiles(signal, control, modeled, raw_peaks, filtered_peaks))

    def scaling(self, mor: MedianOfRatiosNormalization) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
        return mor.scaling_factors({
            "signal": [(self.project, exp.ind) for exp in self.signal],
            "control": [(self.project, exp.ind) for exp in self.control]
        })

    @staticmethod
    def load() -> list["Config"]:
        with open(PKL, 'rb') as stream:
            return pickle.load(stream)
