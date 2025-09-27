from itertools import chain

from biobit.toolkit import nfcore, seqproj

from stories.nextflow.series import internal, SRA

ATTRIBUTE = "title"

folders = sorted(chain(
    SRA.ROOT.glob("PRJ*"),
    internal.ROOT.glob("B*"),
))

for folder in folders:
    print(f"Parsing {folder.name}")
    yaml = folder / "seq-project.yaml"

    if not yaml.exists():
        print("\tNo seq-project.yaml, skipping:", folder.name)
        continue

    project = seqproj.adapter.yaml.load(yaml)

    # Save design as csv
    is_ok = True
    for exp in project.experiments:
        if ATTRIBUTE not in exp.attributes:
            is_ok = False
            print(f"\tSample {exp.sample.ind} is missing {ATTRIBUTE} attribute")

    if is_ok:
        design = folder / "design.csv"
        nfcore.rnaseq.design.from_seqproj(
            project=project,
            saveto=design,
            seqexp2desc=lambda experiment: nfcore.rnaseq.descriptor.from_seqexp(experiment, title_builder=ATTRIBUTE)
        )

    # Generate a new __init__.py file
    with open(folder / "__init__.py", "w") as f:
        f.write("from utils.seqproj import initialize_rnaseq\n\n")
        f.write(f"{folder.name} = initialize_rnaseq(__file__)\n")

    print("Done\n")
