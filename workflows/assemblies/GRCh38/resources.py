from pathlib import Path

ROOT = Path(__file__).parent

RESULTS = ROOT / "results"

FASTA = RESULTS / "GRCh38.primary_assembly.genome.fa.bgz"
GENCODE_GFF3 = RESULTS / "gencode.primary_assembly.annotation.gff3.gz"

GENCODE_VERSION = "47"
