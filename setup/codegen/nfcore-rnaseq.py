import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional, Tuple, Any

from biobit.toolkit import seqproj, nfcore

from resources import (
    HEADER, NFCORE_RNASEQ_DIR, ORGANISM_TO_ASSEMBLY,
    NFCORE_RNASEQ_DESIGN, setup_logging
)

# Template for the auto-generated __init__.py in each project subfolder.
PARSER_FN = """
def {project_name}() -> seqproj.Project:
    from data import {project_name}
    return nfcore.rnaseq.parse.into_seqproj(
        {project_name}(), {results_dir},
        seqexp2descriptor=lambda exp: nfcore.rnaseq.descriptor.from_seqexp(exp, title_builder="title")
    )
"""


def _process_project(folder: Path, data_module: Any) -> Optional[Tuple[str, str, seqproj.Project]]:
    """
    Validates and initializes a single project subdirectory.

    This function handles symlink creation, `__init__.py` generation, and design file
    creation for one project. It returns key project details on success or None on failure.
    """
    logging.info(f"--- Processing project: '{folder.name}' ---")

    # --- 1. Initial Validations ---
    if not hasattr(data_module, folder.name):
        logging.warning(f"Skipping '{folder.name}': No corresponding object in the 'data' module.")
        return None

    init_file = folder / "__init__.py"
    if init_file.exists() and not init_file.read_text().startswith(HEADER):
        logging.error(f"Skipping '{folder.name}': Found a manually-edited '__init__.py'.")
        return None

    try:
        project: seqproj.Project = getattr(data_module, folder.name)()
    except Exception as e:
        logging.error(f"Skipping '{folder.name}': Failed to load its data seq-project. Error: {e}")
        return None

    # --- 2. Determine Genome Assembly ---
    assemblies = {
        ORGANISM_TO_ASSEMBLY[org]
        for smp in project.samples
        for org in smp.organism
        if org in ORGANISM_TO_ASSEMBLY
    }
    if len(assemblies) != 1:
        logging.error(
            f"Skipping '{folder.name}': Expected exactly one assembly, "
            f"but found {len(assemblies)}: {assemblies or 'None'}."
        )
        return None
    assembly = assemblies.pop()
    logging.info(f"Determined assembly for '{folder.name}': {assembly}")

    # --- 3. Filesystem Operations (Symlink and File Generation) ---
    resources_symlink = folder / "resources"
    design_file = folder / NFCORE_RNASEQ_DESIGN

    # Check for a pre-existing directory that would conflict with the symlink
    if resources_symlink.exists() and not resources_symlink.is_symlink():
        logging.error(f"Skipping '{folder.name}': A directory named 'resources' already exists.")
        return None

    try:
        # Create resources symlink
        symlink_target = NFCORE_RNASEQ_DIR.parent / "resources" / "indexes" / assembly
        resources_symlink.unlink(missing_ok=True)  # Remove old symlink if it exists
        relative_target = symlink_target.relative_to(folder, walk_up=True)
        resources_symlink.symlink_to(relative_target, target_is_directory=True)
        logging.info(f"Created symlink: '{resources_symlink}' -> '{relative_target}'")

        # Generate the project's __init__.py
        with open(init_file, "w") as f:
            f.write(HEADER)
            f.write("from pathlib import Path\n\n")
            f.write("from biobit.toolkit import nfcore, seqproj\n\n")
            f.write(f"ROOT = Path(__file__).parent\n\n")
            f.write(PARSER_FN.format(project_name=folder.name, results_dir="ROOT / 'results'"))
        logging.info(f"Generated `__init__.py` for '{folder.name}'.")

        # Create the design file for the nf-core/rnaseq pipeline
        nfcore.rnaseq.design.from_seqproj(
            project=project,
            saveto=design_file,
            seqexp2desc=lambda exp: nfcore.rnaseq.descriptor.from_seqexp(exp, title_builder="title"),
        )
        logging.info(f"Generated design file: '{design_file.name}'")

    except Exception as e:
        logging.error(f"Failed to initialize '{folder.name}' during file generation. Error: {e}")
        # Clean up any partially created files
        resources_symlink.unlink(missing_ok=True)
        init_file.unlink(missing_ok=True)
        design_file.unlink(missing_ok=True)
        return None

    return folder.name, assembly, project


def _generate_toplevel_init(rnaseq_dir: Path, projects_by_assembly: dict):
    logging.info(f"Generating top-level '__init__.py' in '{rnaseq_dir}'...")
    rnaseq_init_file = rnaseq_dir / "__init__.py"

    if rnaseq_init_file.exists() and not rnaseq_init_file.read_text().startswith(HEADER):
        logging.warning(
            "Skipping top-level '__init__.py' generation: an unmanaged file already exists."
        )
        return

    try:
        with open(rnaseq_init_file, "w") as f:
            f.write(HEADER)
            f.write("from pathlib import Path\n\n")

            # Import all projects, sorted alphabetically
            all_projects = sorted(item for projects in projects_by_assembly.values() for item in projects)
            for name, project in all_projects:
                description = f"  # {project.description}" if project.description else ""
                f.write(f"from .{name} import {name}{description}\n")
            f.write("\n")

            f.write("ROOT = Path(__file__).parent\n\n")

            # Create a dictionary grouping project loaders by assembly
            f.write("PER_ASSEMBLY = {\n")
            for assembly, projects in sorted(projects_by_assembly.items()):
                f.write(f"    '{assembly}': [\n")
                for name, _ in sorted(projects):
                    f.write(f"        {name},\n")
                f.write("    ],\n")
            f.write("}\n")
    except IOError as e:
        logging.error(f"Failed to write top-level '__init__.py'. Error: {e}")
        sys.exit(1)


def main():
    if not NFCORE_RNASEQ_DIR.exists():
        logging.error(f"nf-core/rnaseq directory not found: {NFCORE_RNASEQ_DIR}")
        sys.exit(1)

    try:
        import data
    except ImportError as e:
        logging.error(f"Failed to import 'data' module. Ensure DATA_DIR is valid. Error: {e}")
        sys.exit(1)

    by_assembly: dict[str, list[tuple[str, seqproj.Project]]] = defaultdict(list)
    logging.info(f"Scanning subdirectories in '{NFCORE_RNASEQ_DIR}'...")

    # --- Part 1: Process each project subdirectory ---
    for folder in sorted(NFCORE_RNASEQ_DIR.iterdir()):
        if not folder.is_dir():
            continue  # Skip files, symlinks, etc.

        if result := _process_project(folder, data):
            project_name, assembly, project_obj = result
            by_assembly[assembly].append((project_name, project_obj))

    # --- Part 2: Generate the top-level __init__.py ---
    if not by_assembly:
        logging.warning("No valid projects were found to initialize.")
        return

    _generate_toplevel_init(NFCORE_RNASEQ_DIR, by_assembly)
    logging.info("All data modules initialized successfully.")


if __name__ == "__main__":
    setup_logging()
    main()
