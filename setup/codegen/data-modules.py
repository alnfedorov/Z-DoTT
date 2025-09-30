import logging
import sys

from biobit.toolkit import seqproj

from resources import HEADER, DATA_DIR, SEQ_PROJECT_YAML, setup_logging

PARSER_FN = f"""
def {{}}() -> seqproj.Project:
    return seqproj.adapter.yaml.load(ROOT / "{SEQ_PROJECT_YAML}")
"""


def main():
    """
    Scans subdirectories in the DATA_DIR for seq-project.yaml files
    and auto-generates __init__.py modules to make them "importable" from Python.
    """
    if not DATA_DIR.exists():
        logging.error(f"Data directory not found: {DATA_DIR}")
        sys.exit(1)

    projects: dict[str, seqproj.Project] = {}
    logging.info(f"Scanning subdirectories in {DATA_DIR}...")

    # --- Part 1: Generate __init__.py for each project subfolder ---
    for folder in sorted(DATA_DIR.iterdir()):
        if not folder.is_dir():
            continue  # Skip non-directories quietly

        yaml = folder / SEQ_PROJECT_YAML
        if not yaml.exists():
            continue  # Skip folders without the project yaml file

        # Check existing __init__.py before proceeding
        init = folder / "__init__.py"
        if init.exists() and not init.read_text().startswith(HEADER):
            logging.warning(f"Skipping '{folder.name}': found a manually-edited __init__.py.")
            continue

        # Load the project once to avoid redundant I/O
        try:
            project = seqproj.adapter.yaml.load(yaml)
            projects[folder.name] = project
        except Exception as e:
            logging.error(f"Skipping '{folder.name}': failed to parse {SEQ_PROJECT_YAML}. Error: {e}")
            continue

        # Generate the content for the subfolder's __init__.py
        logging.info(f"Generating __init__.py for '{folder.name}'...")
        try:
            with open(init, "w") as f:
                f.write(HEADER)
                f.write("from pathlib import Path\n\n")
                f.write("from biobit.toolkit import seqproj\n\n")
                f.write(f"ROOT = Path(__file__).parent\n\n")
                f.write(PARSER_FN.format(folder.name))
        except Exception as e:
            logging.error(f"Failed to write __init__.py for '{folder.name}'. Error: {e}")
            continue

    # --- Part 2: Generate the top-level data/__init__.py ---
    if not projects:
        logging.warning("No projects found to generate a top-level __init__.py.")
        return

    logging.info(f"Generating top-level __init__.py in {DATA_DIR}...")
    data_init = DATA_DIR / "__init__.py"

    if data_init.exists() and not data_init.read_text().startswith(HEADER):
        logging.warning(
            f"Skipping top-level __init__.py generation: file exists and is not auto-generated."
        )
        return

    try:
        with open(data_init, "w") as f:
            f.write(HEADER)
            f.write("from pathlib import Path\n\n")
            for name, project in sorted(projects.items()):
                description = f"  # {project.description}" if project.description else ""
                f.write(f"from .{name} import {name}{description}\n")
            f.write("\nROOT = Path(__file__).parent\n")
    except Exception as e:
        logging.error(f"Failed to write top-level __init__.py. Error: {e}")
        sys.exit(1)

    logging.info("Data modules initialized successfully.")


if __name__ == "__main__":
    setup_logging()
    main()
