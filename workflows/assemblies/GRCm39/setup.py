import resources
import seqid
from lib import setup, logs
from lib.assembly.annotation.gencode.indexing import index as gencode_index
from lib.assembly.annotation.refseq.indexing import index as refseq_index

logs.setup()

if __name__ == "__main__":
    base = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_{resources.GENCODE_VERSION}"

    # Fasta
    setup.download(
        f"{base}/GRCm39.primary_assembly.genome.fa.gz", resources.FASTA
    )
    setup.rebgzip(resources.FASTA)

    # GENCODE GFF3
    setup.download(
        f"{base}/gencode.v{resources.GENCODE_VERSION}.primary_assembly.annotation.gff3.gz", resources.GENCODE_GFF3
    )
    resources.GENCODE_ANNOTOME.dump(gencode_index("GRCm39", resources.GENCODE_GFF3))

    # RefSeq GFF3 <- Commited directly to the repository
    # setup.download(
    #     "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
    #     "GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz",
    #     resources.REFSEQ_GFF3
    # )
    resources.REFSEQ_ANNOTOME.dump(refseq_index("GRCm39", seqid.from_refseq, resources.REFSEQ_GFF3))
