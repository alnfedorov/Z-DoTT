from pathlib import Path

import pandas as pd
from biobit.toolkit import seqproj

from utils.seqproj import JCCSeq

ROOT = Path(__file__).parent
FASTQ = ROOT / "fastq"


def sample_builder(ind: str, data: pd.DataFrame) -> seqproj.Sample:
    tags = data['Tags'].apply(lambda x: x[1]).unique()
    assert len(tags) == 2 and len({x.replace("input", "").replace("RIP", "") for x in tags}) == 1
    tag = tags[0]

    replica = tag[-1]
    assert replica in {"1", "2", "3", "4"}

    attributes = {"cells": "HT-29", "replica": tag[-1]}

    if tag.startswith("Mock"):
        attributes["treatment"] = "mock"
        organism = {"Homo sapiens"}
    elif tag.startswith("HSV-1"):
        attributes["treatment"] = "HSV-1"
        attributes["HSV-1"] = "F strain[VR-733]"
        attributes["MOI"] = "10"
        attributes["Hours"] = "7"
        organism = {"Homo sapiens", "Herpes simplex virus 1"}
    elif tag.startswith("PR8"):
        attributes["treatment"] = "IAV"
        attributes["IAV"] = "PR8"
        attributes["MOI"] = "5"
        attributes["Hours"] = "7"
        organism = {"Homo sapiens", "Influenza A virus"}
    else:
        raise ValueError(f"Unknown treatment: {tag}")

    return seqproj.Sample(ind=ind, organism=organism, attributes=attributes)


def library_builder(_: str, __: str, data: pd.DataFrame) -> seqproj.Library:
    tags = set(data['Tags'])
    assert len(tags) == 1
    tag = tags.pop()[1][:-2]

    attributes = {}
    if tag.endswith("RIP"):
        tag = tag[:-3]
        if tag.endswith("Z22"):
            rip = "Z22 RIP"
        elif tag.endswith("FLAG"):
            rip = "FLAG RIP"
        else:
            raise ValueError(f"Unknown library: {tags}")
        selection = {"Total RNA", "rRNA depletion", rip}
    elif tag.endswith("input"):
        tag = tag[:-5]
        if tag.endswith("FLAG"):
            attributes["source"] = "FLAG RIP"
        elif tag.endswith("Z22"):
            attributes["source"] = "Z22 RIP"
        else:
            raise ValueError(f"Unknown library: {tags}")
        selection = {"Total RNA", "rRNA depletion"}
    else:
        raise ValueError(f"Unknown library: {tags}")

    return seqproj.Library({"RNA"}, selection, seqproj.Strandedness.Reverse, attributes)


def experiment_builder(
        ind: str, sample: seqproj.Sample, library: seqproj.Library, runs: tuple[seqproj.Run, ...], _: pd.DataFrame
) -> seqproj.Experiment:
    if "FLAG RIP" in library.selection:
        selection = "FLAG-RIP"
    elif "Z22 RIP" in library.selection:
        selection = "Z22-RIP"
    else:
        assert library.selection == {"Total RNA", "rRNA depletion"}, library.selection
        assert "source" in library.attributes, library.attributes
        if library.attributes["source"] == "FLAG RIP":
            selection = "FLAG-input"
        elif library.attributes["source"] == "Z22 RIP":
            selection = "Z22-input"
        else:
            raise ValueError(f"Unknown source: {library.attributes['source']}")

    title = f"{sample.attributes['cells']}_{sample.attributes['treatment']}_{selection}_{sample.attributes['replica']}"
    return seqproj.Experiment(
        ind=ind, sample=sample, library=library, runs=runs, attributes={"title": title}
    )


def project_builder(
        inds: tuple[str, ...], experiments: tuple[seqproj.Experiment, ...], samples: tuple[seqproj.Sample, ...]
):
    assert len(inds) == 1
    return seqproj.Project(f"HT-29 HSV-1/IAV [{ROOT.name}, {inds[0]}]", experiments, samples)


data = JCCSeq.initialize(FASTQ.glob("**/*.fastq.gz"), ROOT)

pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# Remap samples
data["Lane"] = [f"{s}-{l}" for s, l in data[["Sample", "Lane"]].itertuples(index=False)]
sample = data["Tags"].apply(lambda x: x[1])
assert set(sample) == {
    "A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1",
    "A2", "B2", "C2", "D2", "E2", "F2", "G2", "H2",
    "A3", "B3", "C3", "D3", "E3", "F3", "G3", "H3",
    "A4", "B4", "C4", "D4", "E4", "F4", "G4", "H4",
    "A5", "B5", "C5", "D5", "E5", "F5", "G5", "H5",
}

remapping = {
    "A1": "A1+E3", "E3": "A1+E3",
    "B1": "B1+F3", "F3": "B1+F3",
    "C1": "C1+G3", "G3": "C1+G3",
    "D1": "D1+H3", "H3": "D1+H3",
    "E1": "E1+A4", "A4": "E1+A4",
    "F1": "F1+B4", "B4": "F1+B4",
    "G1": "G1+C4", "C4": "G1+C4",
    "H1": "H1+D4", "D4": "H1+D4",
    "A2": "A2+E4", "E4": "A2+E4",
    "B2": "B2+F4", "F4": "B2+F4",
    "C2": "C2+G4", "G4": "C2+G4",
    "D2": "D2+H4", "H4": "D2+H4",
    "E2": "E2+A5", "A5": "E2+A5",
    "F2": "F2+B5", "B5": "F2+B5",
    "G2": "G2+C5", "C5": "G2+C5",
    "H2": "H2+D5", "D5": "H2+D5",
    "A3": "A3+E5", "E5": "A3+E5",
    "B3": "B3+F5", "F5": "B3+F5",
    "C3": "C3+G5", "G5": "C3+G5",
    "D3": "D3+H5", "H5": "D3+H5",
}
data["Sample"] = sample.map(remapping)

assert data['Tags'].apply(lambda x: x[0]).unique() == ["July2024"]
data['Tags'] = data['Tags'].apply(lambda x: x[1:])

project = JCCSeq.parse(
    data,
    seq_machine="Illumina NovaSeq 6000",
    sample_builder=sample_builder,
    library_builder=library_builder,
    experiment_builder=experiment_builder,
    project_builder=project_builder
)
seqproj.adapter.yaml.dump(project, ROOT / "seq-project.yaml")
