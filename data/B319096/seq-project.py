from pathlib import Path

import pandas as pd
from biobit.toolkit import seqproj

from utils.seqproj import JCCSeq

ROOT = Path(__file__).parent
FASTQ = ROOT / "fastq"


def sample_builder(ind: str, data: pd.DataFrame) -> seqproj.Sample:
    tags = set()
    for tag in data["Tags"]:
        assert tag[0] == "Dec2023"
        if tag[-1] == "input":
            tags.add(tag[1:-1])
        else:
            assert tag[-2:] == ("FLAG", "RIP") or tag[-2:] == ("FALG", "RIP")
            tags.add(tag[1:-2])
    assert len(tags) == 1

    treatment, replica = tags.pop()
    attributes = {"cells": "MEF", "replica": replica}

    if treatment == "MCOK":
        treatment = "mock"
        organism = {"Mus musculus"}
    else:
        assert treatment == "HSV-1"
        attributes["HSV-1"] = "F strain[VR-733]"
        organism = {"Mus musculus", "Herpes simplex virus 1"}
    attributes["treatment"] = treatment

    return seqproj.Sample(ind=ind, organism=organism, attributes=attributes)


def library_builder(_: str, __: str, data: pd.DataFrame) -> seqproj.Library:
    tags = set(data['Tags'])
    assert len(tags) == 1
    tags = tags.pop()

    if "input" in tags:
        selection = {"Total RNA", "rRNA depletion"}
    else:
        assert "RIP" in tags and ("FALG" in tags or "FLAG" in tags)
        selection = {"Total RNA", "rRNA depletion", "FLAG RIP"}

    return seqproj.Library({"RNA"}, selection, seqproj.Strandedness.Reverse)


def experiment_builder(
        ind: str, sample: seqproj.Sample, library: seqproj.Library, runs: tuple[seqproj.Run, ...], _: pd.DataFrame
) -> seqproj.Experiment:
    if "FLAG RIP" in library.selection:
        selection = "FLAG-RIP"
    else:
        assert library.selection == {"Total RNA", "rRNA depletion"}
        selection = "input"

    title = f"{sample.attributes['cells']}_{sample.attributes['treatment']}_{selection}_{sample.attributes['replica']}"
    return seqproj.Experiment(
        ind=ind, sample=sample, library=library, runs=runs, attributes={"title": title}
    )


def project_builder(
        inds: tuple[str, ...], experiments: tuple[seqproj.Experiment, ...], samples: tuple[seqproj.Sample, ...]
):
    assert len(inds) == 1
    return seqproj.Project(f"MEF HSV-1 batch 2 [{ROOT.name}, {inds[0]}]", experiments, samples)


data = JCCSeq.initialize(FASTQ.glob("**/*.fastq.gz"), ROOT)
# Fix samples & re-sequenced lanes
remapping = {
    ("S10", "L001"): ("A1+G1", "L001"),
    ("S130", "L001"): ("A1+G1", "L002"),
    ("S16", "L001"): ("A1+G1", "L001"),
    ("S136", "L001"): ("A1+G1", "L002"),
    ("S11", "L001"): ("B1+H1", "L001"),
    ("S131", "L001"): ("B1+H1", "L002"),
    ("S17", "L001"): ("B1+H1", "L001"),
    ("S137", "L001"): ("B1+H1", "L002"),
    ("S12", "L001"): ("C1+A2", "L001"),
    ("S132", "L001"): ("C1+A2", "L002"),
    ("S18", "L001"): ("C1+A2", "L001"),
    ("S138", "L001"): ("C1+A2", "L002"),
    ("S13", "L001"): ("D1+B2", "L001"),
    ("S133", "L001"): ("D1+B2", "L002"),
    ("S19", "L001"): ("D1+B2", "L001"),
    ("S139", "L001"): ("D1+B2", "L002"),
    ("S14", "L001"): ("E1+C2", "L001"),
    ("S134", "L001"): ("E1+C2", "L002"),
    ("S20", "L001"): ("E1+C2", "L001"),
    ("S140", "L001"): ("E1+C2", "L002"),
    ("S15", "L001"): ("F1+D2", "L001"),
    ("S135", "L001"): ("F1+D2", "L002"),
    ("S21", "L001"): ("F1+D2", "L001"),
    ("S141", "L001"): ("F1+D2", "L002")
}
data["Sample"], data["Lane"] = zip(*[remapping[x] for x in data[["Sample", "Lane"]].itertuples(index=False)])

project = JCCSeq.parse(
    data,
    seq_machine="Illumina NovaSeq 6000",
    sample_builder=sample_builder,
    library_builder=library_builder,
    experiment_builder=experiment_builder,
    project_builder=project_builder
)
seqproj.adapter.yaml.dump(project, ROOT / "seq-project.yaml")
