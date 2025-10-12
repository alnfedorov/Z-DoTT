from . import seqid, assembly, resources


def GRCm39() -> assembly.GRCm39:
    return assembly.GRCm39(
        resources.FASTA,
        resources.GENCODE_GFF3, resources.GENCODE_ANNOTOME,
        resources.REFSEQ_GFF3, resources.REFSEQ_ANNOTOME,
        resources.REDIPORTAL,
        resources.REPMASKER, resources.REPMASKER_CLASSIFICATION
    )
