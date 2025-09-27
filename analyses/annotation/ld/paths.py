from pathlib import Path

ROOT = Path(__file__).parent
RESULTS = ROOT / "results"


class gencode:
    root = RESULTS / "gencode"
    pkl = {
        "GRCm39": root / "GRCm39.pkl",
        "CHM13v2": root / "CHM13v2.pkl",
    }


class rna_cores:
    root = RESULTS / "rna-cores"
    pkl = {
        "GRCm39": root / "GRCm39.pkl",
        "CHM13v2": root / "CHM13v2.pkl",
    }


class transcription_boundaries:
    root = RESULTS / "boundaries"
    pkl = {
        "GRCm39": root / "GRCm39.pkl",
        "CHM13v2": root / "CHM13v2.pkl",
    }


class resolved_annotation:
    root = RESULTS / "resolved"
    pkl = {
        "GRCm39": root / "GRCm39.pkl",
        "CHM13v2": root / "CHM13v2.pkl",
    }
