import pickle
from pathlib import Path

PKL = Path(__file__).with_suffix(".pkl")

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import pandas as pd

from utils.repeto import Partition
from stories.RIP import pcalling


@dataclass(frozen=True, slots=True)
class StatTest:
    design: pd.DataFrame  # Design table of each experiment
    formula: str
    alternative: Literal["greaterAbs", "lessAbs", "greater", "less"]
    baselines: dict[str, str]
    comparisons: dict[str, str]

    log2fc: float
    alpha: float

    def __post_init__(self):
        # Sanity checks
        assert not self.design.empty, f"Empty design table: {self}"
        for name, value in self.baselines.items():
            assert value in self.design['group'].values, f"Baseline {value} not found in design table: {self}"

        for name, value in self.comparisons.items():
            column, value = value.split('_', maxsplit=1)
            target, ref = value.split('_vs_')
            assert ref == self.baselines[column], \
                f"Reference group {ref} does not match baseline group {self.baselines[column]}"
            assert target in self.design['group'].values, f"Target value {target} not found in the design table: {self}"
            assert column in self.design.columns, f"Column {column} not found in design table: {self}"


@dataclass(frozen=True, slots=True)
class Config:
    ind: str
    group: str
    elements: tuple[Partition, ...]  # List of groups
    comparisons: tuple[pcalling.Config, ...]  # List of pairwise comparisons
    tests: tuple[StatTest, ...]

    root: Path

    # Annotation results
    alltests: Path = field(init=False)
    sequences: Path = field(init=False)
    localization: Path = field(init=False)
    a2i: Path = field(init=False)
    probes: Path = field(init=False)
    loops: Path = field(init=False)
    summary: Path = field(init=False)

    def __post_init__(self):
        allcmps = set()
        for test in self.tests:
            for cmp in test.comparisons:
                assert cmp not in allcmps, f"Duplicate stat test: {cmp}"
                allcmps.add(cmp)
        object.__setattr__(self, "alltests", self.root / "stat_test.pkl")
        object.__setattr__(self, "summary", self.root / "summary.pkl")
        object.__setattr__(self, "localization", self.root / "localization.pkl")
        object.__setattr__(self, "a2i", self.root / "a2i.pkl")
        object.__setattr__(self, "probes", self.root / "probes.pkl")
        object.__setattr__(self, "loops", self.root / "loops.pkl")
        object.__setattr__(self, "sequences", self.root / "sequences.pkl")

    @staticmethod
    def load() -> list["Config"]:
        with open(PKL, 'rb') as stream:
            return pickle.load(stream)
