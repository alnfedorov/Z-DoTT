from pathlib import Path

from stories.nextflow import series

ROOT = Path(__file__).parent

RESULTS = ROOT.joinpath("results")
VIRAL = ROOT / "viral.bed"


class plots:
    saveto = RESULTS / "plots"
    palette = {
        "Exonic": "#FF7F0E",
        "Intronic": "#1F77B4",
        "Proximal intergenic": "#9467BD",
        "Distal intergenic": "#BD0002"
    }

    mapping = {
        "exon": "Exonic",
        "intron": "Intronic",

        "TES+1kb": "Proximal intergenic",
        "TES+5kb": "Proximal intergenic",
        "TES+10kb": "Proximal intergenic",

        "intergenic": "Distal intergenic",
    }


class reat:
    candidates = RESULTS / "reat" / "candidates"
    filtered = RESULTS / "reat" / "filtered"
    annotated = RESULTS / "reat" / "annotated"
    editome = RESULTS / "reat" / "editome"
    editing_index = RESULTS / "reat" / "editing-index.csv.gz"

    seqproj = RESULTS / "seqproj.pkl"

    experiments = {
        "CHM13v2": series.internal.by_assembly.get("CHM13v2", []) + [
            series.SRA.PRJNA256013,
            series.SRA.PRJNA382632,
        ],
        "GRCm39": series.internal.by_assembly.get("GRCm39", []) + series.SRA.by_assembly.get("GRCm39", [])
    }


class tracks:
    all_annotated = RESULTS / "esites" / "all-annotated"
    all_passed = RESULTS / "esites" / "all-passed"
    expwise = RESULTS / "esites" / "expwise"
    editome = RESULTS / "esites" / "editome"


class thresholds:
    min_coverage = 8
    freqthr = (0.05, 0.45)
    min_edits = 1
    min_samples = 3
