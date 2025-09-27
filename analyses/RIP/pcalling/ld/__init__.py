import gzip
from pathlib import Path

from .config import Config

ROOT = Path(__file__).parent
RESULTS = ROOT / "results"

EFFECTIVE_GENOME_SIZE = RESULTS / "effective-genome-size.pkl"
RNA_MODELS = RESULTS / "rna-models.pkl"


def write_gz(path: str, lines: list[str]):
    with gzip.open(path, 'wt') as stream:
        stream.writelines(lines)
