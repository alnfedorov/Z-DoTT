from pathlib import Path

from biobit.toolkit import seqproj, nfcore

FOLDER = Path(__file__).parent

FASTQ = FOLDER / "fastq"

loaded = nfcore.fetchngs.load_seqproj(
    samplesheet=FASTQ / "samplesheet.csv",
    fastq_root=Path("fastq")
)
assert len(loaded) == 1, f"Expected 1 seqproj, got {len(loaded)} for {FOLDER.name}"
project = loaded[0]

titles = set()
for exp in project.experiments:
    assert exp.sample.attributes["title"] == exp.sample.description
    title = exp.sample.attributes.pop("title")

    donor, virus, moi, timepoint, replica = title.split("_")
    timepoint = timepoint.lstrip("0")
    donor = donor.replace("MDM-", "")
    treatment = {"MOI2": virus, "MOI0": "mock"}[moi]

    object.__setattr__(exp.sample, "organism", {"Homo sapiens", f"IAV ({virus})"})
    object.__setattr__(exp.library, "strandedness", seqproj.Strandedness.Reverse)

    exp.sample.attributes["cells"] = "MDM"
    exp.sample.attributes["donor"] = donor
    exp.sample.attributes["IAV"] = treatment
    exp.sample.attributes["time-point"] = timepoint
    exp.sample.attributes["replica"] = replica

    title = "_".join([
        "MDM", exp.sample.attributes["donor"], exp.sample.attributes["IAV"],
        exp.sample.attributes["time-point"], exp.sample.attributes["replica"]
    ])
    assert title not in titles, title
    exp.attributes["title"] = title
    titles.add(title)

seqproj.adapter.yaml.dump(project, FOLDER / "seq-project.yaml")
