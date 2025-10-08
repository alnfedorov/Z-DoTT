from . import assembly, resources
from .. import HSV1, CHM13v2, IAV


def zdhs() -> assembly.zdhs:
    return assembly.zdhs(resources.FASTA, resources.GFF3, HSV1(), IAV(), CHM13v2())
