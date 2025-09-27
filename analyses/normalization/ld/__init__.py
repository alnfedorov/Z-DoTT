from pathlib import Path

from stories import nextflow
from .normbin import NormBin

ROOT = Path(__file__).parent
RESULTS = ROOT / "results"

# Counting results
BINS = RESULTS / "bins.pkl"
METRICS = RESULTS / "metrics.pkl"
FRAGMENTS = RESULTS / "fragments-distribution.pkl"
VIRAL_RNA_LOAD = RESULTS / "viral-rna-load.pkl"

SERIES = {
    "GRCm39": nextflow.series.internal.by_assembly["GRCm39"],
    "CHM13v2": [
        *nextflow.series.internal.by_assembly["CHM13v2"],
        *nextflow.series.SRA.by_assembly["CHM13v2"]
    ]
}
