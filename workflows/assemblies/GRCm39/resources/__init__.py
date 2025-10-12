from pathlib import Path

from lib.pkl import PklData

ROOT = Path(__file__).parent

REDIPORTAL = ROOT / "rediportal.GRCm39.bed.gz"
REPMASKER = ROOT / "repmasker.GRCm39.bed.gz"
REPMASKER_CLASSIFICATION = ROOT / "repmasker.GRCm39.classification.tsv.gz"

RESULTS = ROOT.parent / "results"

FASTA = RESULTS / "GRCm39.primary_assembly.genome.fa.bgz"
GENCODE_GFF3 = RESULTS / "gencode.primary_assembly.annotation.gff3.gz"
GENCODE_ANNOTOME = PklData(RESULTS / "gencode.primary_assembly.annotation.pkl")
REFSEQ_GFF3 = RESULTS / "GCF_000001635.27_GRCm39_genomic.gff3.gz"
REFSEQ_ANNOTOME = PklData(RESULTS / "GCF_000001635.27_GRCm39_genomic.pkl")

GENCODE_VERSION = "M36"
