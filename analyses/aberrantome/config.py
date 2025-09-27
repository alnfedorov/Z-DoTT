from dataclasses import dataclass, field
from pathlib import Path

from biobit.toolkit.seqproj import Experiment


@dataclass(frozen=True)
class Config:
    ind: str
    assembly: str

    # Project ind + Experiment
    treatment: tuple[tuple[str, Experiment], ...]
    control: tuple[tuple[str, Experiment], ...]

    # Labels
    treatment_label: str
    control_label: str

    # Path to a results folder
    results: Path

    # Path to estimated scores and results of statistical tests (csv tables)
    read_through: Path = field(init=False)
    read_in: Path = field(init=False)
    divergent: Path = field(init=False)
    intronic: Path = field(init=False)

    def __post_init__(self):
        # Ensure that the treatment and control are sorted
        treatment = tuple(sorted(self.treatment, key=lambda x: (x[0], x[1].ind)))
        control = tuple(sorted(self.control, key=lambda x: (x[0], x[1].ind)))

        uniqtrt = set((x[0], x[1].ind) for x in treatment)
        uniqctrl = set((x[0], x[1].ind) for x in control)
        assert len(uniqtrt & uniqctrl) == 0, f"Treatment and control must be disjoint, got {uniqtrt & uniqctrl}"

        results = self.results / self.ind
        read_through = results / "read_through.csv.gz"
        read_in = results / "read_in.csv.gz"
        divergent = results / "divergent.csv.gz"
        intronic = results / "intronic.csv.gz"

        for obj, ket in [
            (treatment, "treatment"), (control, "control"), (results, "results"), (divergent, "divergent"),
            (read_through, "read_through"), (read_in, "read_in"), (intronic, "intronic")
        ]:
            object.__setattr__(self, ket, obj)

    def treatment_key(self) -> tuple[tuple[str, str], ...]:
        return tuple(sorted((x[0], x[1].ind) for x in self.treatment))

    def control_key(self) -> tuple[tuple[str, str], ...]:
        return tuple(sorted((x[0], x[1].ind) for x in self.control))
