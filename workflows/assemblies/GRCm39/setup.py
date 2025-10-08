import resources
from lib import setup

if __name__ == "__main__":
    base = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_{resources.GENCODE_VERSION}"

    # Fasta
    setup.download(
        f"{base}/GRCm39.primary_assembly.genome.fa.gz", resources.FASTA
    )
    setup.ungzip_then_bgzip(resources.FASTA)

    # GENCODE GFF3
    setup.download(
        f"{base}/gencode.v{resources.GENCODE_VERSION}.primary_assembly.annotation.gff3.gz", resources.GENCODE_GFF3
    )

    # RefSeq GFF3
    setup.download(
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
        "GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz",
        resources.REFSEQ_GFF3
    )
