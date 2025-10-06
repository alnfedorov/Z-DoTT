import logging
from pathlib import Path

HEADER = "# AUTO-GENERATED FILE - DO NOT EDIT\n"


SEQ_PROJECT_YAML = "seq-project.yaml"
NFCORE_RNASEQ_DESIGN = "design.csv"

ORGANISM_TO_ASSEMBLY = {
    "Homo sapiens": "CHM13v2", "Mus musculus": "GRCm39"
}


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
