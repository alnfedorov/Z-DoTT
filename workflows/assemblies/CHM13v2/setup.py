import logging

import lib.assembly.annotation.gencode.liftoff.indexing
import lib.assembly.annotation.refseq.indexing
import resources
import seqid
from lib import logs, setup

logs.setup()

if __name__ == "__main__":
    # Fasta
    setup.download(
        "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
        resources.FASTA
    )
    setup.rebgzip(resources.FASTA)

    # # RefSeq GFF3 <- Commited to the repository
    # setup.download(
    #     "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/"
    #     "GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz",
    #     resources.REFSEQ_GFF3
    # )
    resources.REFSEQ_ANNOTOME.dump(
        lib.assembly.annotation.refseq.indexing.index(
            "CHM13v2", seqid.from_refseq, resources.REFSEQ_GFF3
        )
    )

    if resources.GENCODE_LIFTOFF_GFF3.exists():
        resources.GENCODE_LIFTOFF_ANNOTOME.dump(
            lib.assembly.annotation.gencode.liftoff.indexing.index("CHM13v2", "GRCh38", resources.GENCODE_LIFTOFF_GFF3)
        )
    else:
        logging.warning(
            "To complete the setup, liftoff the GENCODE annotation from GRCh38 to CHM13v2 and re-run setup.py"
        )
