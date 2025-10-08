import logging

import resources
from lib import setup

if __name__ == "__main__":
    # Fasta
    setup.download(
        "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
        resources.FASTA
    )
    setup.ungzip_then_bgzip(resources.FASTA)

    # RefSeq GFF3
    setup.download(
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/"
        "GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz",
        resources.REFSEQ_GFF3
    )

    logging.warning("To complete the setup, liftoff the GENCODE annotation from GRCh38 to CHM13v2.0")
