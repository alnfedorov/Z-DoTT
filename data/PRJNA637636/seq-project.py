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

for exp in project.experiments:
    assert exp.sample.attributes["title"] == exp.sample.description
    title = exp.sample.attributes.pop("title")

    if title.endswith("Rep. 1"):
        exp.sample.attributes["replica"] = "1"
    else:
        assert title.endswith("Rep. 2")
        exp.sample.attributes["replica"] = "2"

    title = title[:-len("Rep. 1")].strip()
    virus = {
        "Mock": "mock",
        "WT strain 17": "Strain 17",
        "delta vhs": "delta VHS",
        "TsK": "TsK",
        "delta ICP27": "delta ICP27",
        "delta ICP4": "delta ICP4",
        "delta ICP22": "delta ICP22",
        "delta IPC0": "delta ICP0",
    }[title]
    exp.sample.attributes["HSV-1"] = virus

    time_window = {
        "Mock": "mock",
        "WT strain 17": "7-8hpi",
        "delta vhs": "7-8hpi",
        "TsK": "7-8hpi",
        "delta ICP27": "7-8hpi",
        "delta ICP4": "7-8hpi",
        "delta ICP22": "11-12hpi",
        "delta IPC0": "11-12hpi",
    }[title]
    exp.sample.attributes["time-point"] = time_window

    exp.sample.attributes["cells"] = "HFF"
    exp.sample.attributes["protocol"] = "4sU-RNA"

    object.__setattr__(exp.sample, "organism", {"Homo sapiens", "HSV-1"})
    object.__setattr__(exp.library, "selection", {"4sU", "rRNA depletion"})
    object.__setattr__(exp.library, "strandedness", seqproj.Strandedness.Reverse)

    title = "_".join([
        exp.sample.attributes["cells"], exp.sample.attributes["HSV-1"].replace(" ", "-"),
        exp.sample.attributes["time-point"], exp.sample.attributes["replica"]
    ]).replace("mock_mock", "mock")
    exp.attributes["title"] = title
    print(title)

seqproj.adapter.yaml.dump(project, FOLDER / "seq-project.yaml")
