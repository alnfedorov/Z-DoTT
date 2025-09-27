from pathlib import Path

from stories import nextflow

ROOT = Path(__file__).parent

RESULTS = ROOT / "results"

PCA = RESULTS / "pca"

ANNOTATION = {
    "GRCm39": RESULTS / "annotation" / "GRCm39.pkl",
    "GRCh38": RESULTS / "annotation" / "GRCh38.pkl",
    "CHM13v2": RESULTS / "annotation" / "CHM13v2.pkl",
}


class BIOTYPES:
    ROOT = RESULTS / "biotypes"

    PLOTS = ROOT / "plots"

    METRICS = ROOT / "metrics.pkl"
    RAW_FRAGMENTS = ROOT / "raw-fragments-distribution.pkl"
    COLLAPSED_FRAGMENTS = ROOT / "collapsed-fragments-distribution.pkl"

    series = {
        "GRCm39": [
            nextflow.series.internal.B256178,  # IAV infection of MEFs (batch 1)
            nextflow.series.internal.B261790,  # HSV-1 infection of MEFs (batch 1)
            nextflow.series.internal.B319096,  # HSV-1 infection of MEFs (batch 2)
        ],
        "CHM13v2": [
            nextflow.series.internal.B831009,  # HSV-1/IAV RIP from HT-29s
        ]
    }
