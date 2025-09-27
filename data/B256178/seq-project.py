from pathlib import Path

import pandas as pd
from biobit.toolkit import seqproj

from utils.seqproj import JCCSeq

ROOT = Path(__file__).parent
FASTQ = ROOT / "fastq"


def sample_builder(ind: str, data: pd.DataFrame) -> seqproj.Sample:
    tags = set(data["Tags"])
    assert len(tags) == 2 and all(x[0] == "Feb2022" and x[1] == "IAV" for x in tags), tags
    subind1, subind2 = tags.pop()[-1], tags.pop()[-1]
    assert ind == f"{subind1}+{subind2}" or ind == f"{subind2}+{subind1}", (ind, subind1, subind2)

    attributes = {"cells": "MEF"}
    if {subind1, subind2}.issubset({"A1", "B1", "C1", "G1", "H1", "A2", "E2", "F2", "G2", "C3", "D3", "E3"}):
        attributes["treatment"] = "mock"
        organism = {"Mus musculus"}
    elif {subind1, subind2}.issubset({"D1", "E1", "F1", "B2", "C2", "D2", "H2", "A3", "B3", "F3", "G3", "H3"}):
        attributes["treatment"] = "IAV"
        attributes["time-point"] = "8h"
        attributes["IAV"] = "PR8"
        organism = {"Mus musculus", "Influenza A virus"}
    else:
        raise ValueError(f"Unknown plates: {(ind, subind1, subind2)}")

    attributes["replica"] = {
        "A1+E2": "1", "B1+F2": "1", "C1+G2": "1", "D1+H2": "1", "E1+A3": "1", "F1+B3": "1", "G1+C3": "2", "H1+D3": "2",
        "A2+E3": "2", "B2+F3": "2", "C2+G3": "2", "D2+H3": "2",
    }[ind]

    return seqproj.Sample(ind=ind, organism=organism, attributes=attributes)


def library_builder(_: str, __: str, data: pd.DataFrame) -> seqproj.Library:
    tags = set(data['Tags'])
    assert len(tags) == 1
    _, _, plate = tags.pop()

    if plate in {"E2", "F2", "G2", "H2", "A3", "B3", "C3", "D3", "E3", "F3", "G3", "H3"}:
        # Input tags
        selection = {"Total RNA", "rRNA depletion"}
    elif plate in {"A1", "D1", "G1", "B2"}:
        # IgG tags
        selection = {"Total RNA", "rRNA depletion", "IgG RIP"}
    elif plate in {"B1", "E1", "H1", "C2"}:
        selection = {"Total RNA", "rRNA depletion", "Z22 RIP"}
    elif plate in {"C1", "F1", "A2", "D2"}:
        selection = {"Total RNA", "rRNA depletion", "FLAG RIP"}
    else:
        raise ValueError(f"Unknown plate: {plate}")

    return seqproj.Library({"RNA"}, selection, seqproj.Strandedness.Reverse)


def experiment_builder(
        ind: str, sample: seqproj.Sample, library: seqproj.Library, runs: tuple[seqproj.Run, ...], _: pd.DataFrame
) -> seqproj.Experiment:
    if "FLAG RIP" in library.selection:
        selection = "FLAG-RIP"
    elif "Z22 RIP" in library.selection:
        selection = "Z22-RIP"
    elif "IgG RIP" in library.selection:
        selection = "IgG-RIP"
    else:
        assert library.selection == {"Total RNA", "rRNA depletion"}
        selection = "input"

    title = f"{sample.attributes['cells']}_{sample.attributes['treatment']}_{selection}_{sample.attributes['replica']}"
    return seqproj.Experiment(ind=ind, sample=sample, library=library, runs=runs, attributes={"title": title})


def project_builder(
        inds: tuple[str, ...], experiments: tuple[seqproj.Experiment, ...], samples: tuple[seqproj.Sample, ...]
):
    assert len(inds) == 1
    return seqproj.Project(f"MEF IAV batch 1 [{ROOT.name}, {inds[0]}]", experiments, samples)


data = JCCSeq.initialize(FASTQ.glob("**/*.fastq.gz"), ROOT)

# Remap input & RIP samples
remapping = {
    'A1': 'A1+E2', 'B1': 'B1+F2', 'C1': 'C1+G2', 'D1': 'D1+H2', 'E1': 'E1+A3', 'F1': 'F1+B3', 'G1': 'G1+C3',
    'H1': 'H1+D3', 'A2': 'A2+E3', 'B2': 'B2+F3', 'C2': 'C2+G3', 'D2': 'D2+H3', 'E2': 'A1+E2', 'F2': 'B1+F2',
    'G2': 'C1+G2', 'H2': 'D1+H2', 'A3': 'E1+A3', 'B3': 'F1+B3', 'C3': 'G1+C3', 'D3': 'H1+D3', 'E3': 'A2+E3',
    'F3': 'B2+F3', 'G3': 'C2+G3', 'H3': 'D2+H3',
}
data["Lane"] = [f"{s}-{l}" for s, l in data[["Sample", "Lane"]].itertuples(index=False)]
data["Sample"] = data["Tags"].apply(lambda x: x[2]).map(remapping)
assert data.isnull().sum().sum() == 0, data.isnull().sum()

project = JCCSeq.parse(
    data,
    seq_machine="Illumina NovaSeq 6000",
    sample_builder=sample_builder,
    library_builder=library_builder,
    experiment_builder=experiment_builder,
    project_builder=project_builder
)
seqproj.adapter.yaml.dump(project, ROOT / "seq-project.yaml")
