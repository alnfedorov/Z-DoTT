from pathlib import Path

ROOT = Path(__file__).parent

FASTA = ROOT / "sequence.fa"
CUSTOM_ANNOTATION_GFF3 = ROOT / "sequence.gff3"  # doi.org-10.1038-s41467-020-15992-5 + NCBI LAT RNA + manual fixes
NCBI_ANNOTATION_GFF3 = ROOT / "annotation.ncbi.gff3"
