import logging
from pathlib import Path

from biobit.toolkit import seqproj

logging.basicConfig(level=logging.INFO)

DATA = Path(__file__).parent.parent.parent / "data"
assert DATA.exists(), f"Data directory not found: {DATA}"

SEQ_PROJECT_YAML = "seq-project.yaml"
HEADER = "# AUTO-GENERATED FILE - DO NOT EDIT\n"

projects: dict[str, seqproj.Project] = {}
for folder in DATA.iterdir():
    logging.info(f"Checking {folder}")
    if not folder.is_dir():
        logging.info(f"\tSkipping {folder} - not a directory")
        continue

    yaml = folder / SEQ_PROJECT_YAML
    if not yaml.exists():
        logging.info(f"\tSkipping {folder} - no seq-project.yaml found")
        continue

    logging.info(f"Processing {folder}")
    init = folder / "__init__.py"
    if init.exists():
        content = init.read_text()
        if not content.startswith(HEADER):
            logging.info(f"\tSkipping {folder} - __init__.py exists and is not auto-generated")
            continue

    with open(init, "w") as f:
        f.write(HEADER)
        f.write("from pathlib import Path\n\n")
        f.write("from biobit.toolkit import seqproj\n\n")
        f.write("ROOT = Path(__file__).parent\n")
        f.write(f"{folder.name} = seqproj.adapter.yaml.load(ROOT / \"{yaml.name}\")\n")
    logging.info("\tWritten __init__.py")

    projects[folder.name] = seqproj.adapter.yaml.load(yaml)

logging.info("Generating data/__init__.py")
init = DATA / "__init__.py"

process = True
if init.exists():
    content = init.read_text()
    if not content.startswith(HEADER):
        logging.info(f"\tSkipping data/__init__.py - exists and is not auto-generated")
        process = False
if process:
    with open(init, "w") as f:
        f.write(HEADER)
        f.write("from pathlib import Path\n\n")
        for prjind in sorted(projects.keys()):
            description = projects[prjind].description
            description = f" # {description}" if description else ""
            f.write(f"from .{prjind} import {prjind}{description}\n")
        f.write("\n")
        f.write("ROOT = Path(__file__).parent\n")
