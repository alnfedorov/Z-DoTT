import resources
from lib import setup

if __name__ == "__main__":
    base = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{resources.GENCODE_VERSION}"
    setup.download(
        f"{base}/GRCh38.primary_assembly.genome.fa.gz",
        resources.FASTA
    )
    setup.ungzip_then_bgzip(resources.FASTA)
    setup.index_fasta(resources.FASTA)

    setup.download(
        f"{base}/gencode.v{resources.GENCODE_VERSION}.primary_assembly.annotation.gff3.gz",
        resources.GENCODE_GFF3
    )
