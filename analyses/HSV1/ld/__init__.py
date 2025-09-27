from pathlib import Path

from biobit.core.loc import PerStrand

from . import visual, features, enrichment

ROOT = Path(__file__).parent
RESULTS = ROOT / "results"

PLOTS = RESULTS / "plots"
ENRICHMENT = RESULTS / "enrichment"

MIN_ENRICHMENT = 2
MAX_ENRICHMENT_DIFF = 0.25
BINSIZE = 128


class annotation:
    root = RESULTS / "annotation"

    primers = root / "primers.bed"
    miRNA = root / "miRNA.bed"
    repeats = root / "repeats.bed"
    orfs = PerStrand(forward=root / "orfs.fwd.bed", reverse=root / "orfs.fwd.bed")
    genes = PerStrand(forward=root / "genes.fwd.bed", reverse=root / "genes.rev.bed")

    circos = root / "circos.pkl"


class palette:
    forward = "#f0615b"
    reverse = "#617fbf"
    terminal = "#666666"


class MFE:
    root = RESULTS / "MFE"

    # All sequences with ZH score within top X% are considered Z-prone
    zh_quantile = 0.95
    # Min/max stem length in dinucleotides. Stems longer than maxdn are windowed as maxdn sub-stems
    dnlim = (3, 6)
