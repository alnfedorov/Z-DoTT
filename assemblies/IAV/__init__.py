from pathlib import Path

ROOT = Path(__file__).parent

# Note: The IAV sequence was corrected manually for observed variations
fasta = ROOT / "sequence.fa"
gff3 = ROOT / "sequence.gff3"

name = "IAV"
organism = "Influenza A virus"

segments = {
    "NC_002023.1": 2341,
    "NC_002022.1": 2233,
    "NC_002021.1": 2341,
    "NC_002020.1": 890,
    "NC_002019.1": 1565,
    "NC_002018.1": 1413,
    "NC_002017.1": 1778,
    "NC_002016.1": 1027,
}

__all__ = ["fasta", "gff3", "segments"]
