from . import resources, assembly


def GRCh38() -> assembly.GRCh38:
    return assembly.GRCh38(resources.FASTA, resources.GENCODE_GFF3)
