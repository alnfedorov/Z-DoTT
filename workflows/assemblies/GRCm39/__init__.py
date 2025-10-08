from . import seqid, assembly, resources


def GRCm39() -> assembly.GRCm39:
    return assembly.GRCm39(
        resources.FASTA, resources.GENCODE_GFF3, resources.REFSEQ_GFF3, resources.REDIPORTAL, resources.REPMASKER,
        resources.REPMASKER_CLASSIFICATION
    )
