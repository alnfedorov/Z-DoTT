from pathlib import Path

from . import utils

ROOT = Path(__file__).parent

fasta = ROOT / "sequence.fa"
# doi.org-10.1038-s41467-020-15992-5 + NCBI LAT RNA + manual fixes
gff3 = ROOT / "sequence.gff3"
ncbi = ROOT / "annotation.ncbi.gff3"

name = "HSV-1"
organism = "Herpes simplex virus 1"

contig = "NC_001806.2"
size = 152222
segments = {"NC_001806.2": 152222}

nagnag_splicing = {"TRL2-iso2", "UL36.6-iso2", "IRL2-iso2"}
terminal_repeats = [
    (0, 9213), (145589, 152222),
]

__all__ = ["utils", "fasta", "gff3", "ncbi", "contig", "size", "segments", "nagnag_splicing", "terminal_repeats"]
