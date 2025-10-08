from . import seqid, assembly, resources


def CHM13v2() -> assembly.CHM13v2:
    return assembly.CHM13v2(
        resources.FASTA, resources.GENCODE_LIFTOFF_GFF3, resources.REFSEQ_GFF3, resources.REDIPORTAL,
        resources.REPMASKER, resources.REPMASKER_CLASSIFICATION
    )
