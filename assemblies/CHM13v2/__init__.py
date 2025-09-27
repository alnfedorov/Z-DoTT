from pathlib import Path

from biobit.deprecated.repmasker import RepmaskerClassification

from . import seqid, refseq, gencode

ROOT = Path(__file__).parent

name = "CHM13v2"
organism = "Homo sapiens"

# Sequence
fasta = ROOT / "CHM13v2.fa.bgz"

# RepeatMasker
repcls = RepmaskerClassification(ROOT / "repmasker.CHM13v2.classification.tsv.gz")
repmasker = ROOT / "repmasker.CHM13v2.bed.gz"

# REDI portal
rediportal = ROOT / "rediportal.CHM13v2.bed.gz"
