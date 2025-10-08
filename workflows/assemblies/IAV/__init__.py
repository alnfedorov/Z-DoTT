from . import resources, assembly


def IAV() -> assembly.IAV:
    # Note: The IAV sequence was corrected for observed variations
    return assembly.IAV(resources.FASTA, resources.ANNOTATION_GFF3)
