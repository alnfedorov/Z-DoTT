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
    title = title.replace(" hpi", "hpi")

    object.__setattr__(exp.sample, "organism", {"Homo sapiens", "HSV-1 (strain 17)"})
    exp.sample.attributes["cells"] = "HFF"

    if "4sU-RNA" in title:
        object.__setattr__(exp.library, "selection", {"4sU", })

        title = title.replace("4sU-RNA - ", "")
        timewindow, rep, repn = title.split()
        exp.sample.attributes["4sU"] = "500 uM"
        protocol = "4sU-RNA"
    else:
        assert "Total RNA" in title
        object.__setattr__(exp.library, "selection", {"Total RNA", "rRNA depletion"})
        title = title.replace("Total RNA - ", "")
        timewindow, rep, repn = title.split()
        protocol = "Total-RNA"

    repn = int(repn)
    assert rep == "Rep"
    assert repn in (1, 2, 3)
    exp.sample.attributes["time-point"] = timewindow
    exp.sample.attributes["protocol"] = protocol
    exp.sample.attributes["replica"] = str(repn)

    object.__setattr__(exp.library, "strandedness", seqproj.Strandedness.Reverse)

    title = "_".join([
        exp.sample.attributes["cells"], exp.sample.attributes["time-point"], protocol, exp.sample.attributes["replica"]
    ])
    exp.attributes["title"] = title

seqproj.adapter.yaml.dump(project, FOLDER / "seq-project.yaml")
