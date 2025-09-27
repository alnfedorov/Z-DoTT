from pathlib import Path

from biobit.toolkit import seqproj, nfcore

FOLDER = Path(__file__).parent

FASTQ = FOLDER / "fastq"

loaded = nfcore.fetchngs.load_seqproj(
    samplesheet=FASTQ / "samplesheet.csv",
    fastq_root=Path("fastq")
)
assert len(loaded) == 1, f"Expected 1 seqproject, got {len(loaded)} for {FOLDER.name}"
project = loaded[0]

forward_layout = {
    "ERX12477976",
    "ERX12477981",
    "ERX12477982",
    "ERX12477983",
    "ERX12477985",
    "ERX12477987",
    "ERX12477989",
    "ERX12477990",
    "ERX12477994",
    "ERX12478004",
    "ERX12478005",
    "ERX12478009",
}

for exp in project.experiments:
    assert exp.sample.attributes["title"] == exp.sample.description
    title = exp.sample.attributes.pop("title")

    hpi, condition, fraction, replica = title.split("_")
    condition = {
        "mock": "mock",
        "IAV": "IAV-WT",
        "dNS1": "IAV-dNS1",
    }[condition]
    assert hpi in {"8hpi", "16hpi"}
    assert fraction in {"nucleus", "cytosol"}
    assert replica in {"1", "2", "3"}

    exp.sample.attributes['cells'] = "A549"
    if condition != "mock":
        exp.sample.attributes['strain'] = "A/WSN/1933"
    exp.sample.attributes['IAV'] = condition
    exp.sample.attributes['time-point'] = hpi
    exp.sample.attributes['fraction'] = fraction
    exp.sample.attributes['replica'] = replica

    object.__setattr__(exp.sample, "organism", {"Homo sapiens", "IAV"})
    object.__setattr__(exp.library, "selection", {"rRNA depletion"})

    if exp.ind in forward_layout:
        object.__setattr__(exp.library, "strandedness", seqproj.Strandedness.Forward)
        forward_layout.remove(exp.ind)
    else:
        object.__setattr__(exp.library, "strandedness", seqproj.Strandedness.Reverse)

    # A549_IAV-WT_

    title = "_".join([
        exp.sample.attributes[k] for k in ("cells", "IAV", "time-point", "fraction", "replica")
    ])
    exp.attributes["title"] = title
    print(title)
assert len(forward_layout) == 0, forward_layout

seqproj.adapter.yaml.dump(project, FOLDER / "seq-project.yaml")

