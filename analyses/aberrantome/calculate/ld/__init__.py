import pickle
from pathlib import Path

from . import thresholds as thr

ROOT = Path(__file__).parent

RESULTS = ROOT / "results"

# RNA models for each assembly
RNA = RESULTS / "rna.pkl"


class counts:
    root = RESULTS / "counts"

    metrics = root / "metrics.pkl"

    read_through = root / "read-through.pkl"
    read_in = root / "read-in.pkl"
    divergent = root / "divergent.pkl"
    intronic = root / "intronic.pkl"


class comparisons:
    from stories import nextflow

    series = {
        "CHM13v2": [
            nextflow.series.internal.B831009,
            nextflow.series.SRA.PRJEB75711,
            nextflow.series.SRA.PRJNA256013,
            nextflow.series.SRA.PRJNA637636,
        ],
        "GRCm39": [
            nextflow.series.internal.B256178,
            nextflow.series.internal.B261790
        ]
    }

    pkl = RESULTS / "comparisons.pkl"

    @staticmethod
    def load():
        with open(comparisons.pkl, "rb") as f:
            return pickle.load(f)
