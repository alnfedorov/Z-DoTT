from pathlib import Path

ROOT = Path(__file__).parent

REDIPORTAL = ROOT / "rediportal.GRCm39.bed.gz"
REPMASKER = ROOT / "repmasker.GRCm39.bed.gz"
REPMASKER_CLASSIFICATION = ROOT / "repmasker.GRCm39.classification.tsv.gz"

RESULTS = ROOT.parent / "results"

FASTA = RESULTS / "GRCm39.primary_assembly.genome.fa.bgz"
GENCODE_GFF3 = RESULTS / "gencode.primary_assembly.annotation.gff3.gz"
REFSEQ_GFF3 = RESULTS / "GCF_000001635.27_GRCm39_genomic.gff3.gz"

GENCODE_VERSION = "M36"
