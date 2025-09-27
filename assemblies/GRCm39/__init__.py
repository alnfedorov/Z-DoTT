from pathlib import Path

from biobit.deprecated.repmasker import RepmaskerClassification

from . import seqid, gencode, refseq

ROOT = Path(__file__).parent

name = "GRCm39"
organism = "Mus musculus"

# Sequence
fasta = ROOT / "GRCm39.primary_assembly.genome.fa.bgz"

# RepeatMasker
repcls = RepmaskerClassification(ROOT / "repmasker.GRCm39.classification.tsv.gz")
repmasker = ROOT / "repmasker.GRCm39.bed.gz"

# REDI portal
rediportal = ROOT / "rediportal.GRCm39.bed.gz"
