from . import assembly, resources


def HSV1() -> assembly.HSV1:
    return assembly.HSV1(resources.FASTA, resources.CUSTOM_ANNOTATION_GFF3)
