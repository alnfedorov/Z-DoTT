from pathlib import Path

from lib.pkl import PklData

ROOT = Path(__file__).parent

REDIPORTAL = ROOT / "rediportal.CHM13v2.bed.gz"
REPMASKER = ROOT / "repmasker.CHM13v2.bed.gz"
REPMASKER_CLASSIFICATION = ROOT / "repmasker.CHM13v2.classification.tsv.gz"

RESULTS = ROOT.parent / "results"

FASTA = RESULTS / "CHM13v2.fa.bgz"
REFSEQ_GFF3 = RESULTS / "CHM13v2.refseq.gff.gz"
REFSEQ_ANNOTOME = PklData(RESULTS / "CHM13v2.refseq.pkl")
GENCODE_LIFTOFF_GFF3 = RESULTS / "CHM13v2.liftoff+gencode.gff3.gz"
GENCODE_LIFTOFF_ANNOTOME = PklData(RESULTS / "CHM13v2.liftoff+gencode.pkl")
