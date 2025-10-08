from pathlib import Path

ROOT = Path(__file__).parent

REDIPORTAL = ROOT / "rediportal.CHM13v2.bed.gz"
REPMASKER = ROOT / "repmasker.CHM13v2.bed.gz"
REPMASKER_CLASSIFICATION = ROOT / "repmasker.CHM13v2.classification.tsv.gz"

RESULTS = ROOT.parent / "results"

FASTA = RESULTS / "CHM13v2.fa.bgz"
REFSEQ_GFF3 = RESULTS / "CHM13v2.refseq.gff.gz"
GENCODE_LIFTOFF_GFF3 = RESULTS / "CHM13v2.liftoff+gencode.gff3.gz"
