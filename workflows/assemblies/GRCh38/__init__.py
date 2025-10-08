import assembly
import resources


def GRCh38() -> assembly.GRCh38:
    return assembly.GRCh38(resources.FASTA, resources.GENCODE_GFF3)
