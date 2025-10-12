from . import assembly, resources
from .. import HSV1, GRCm39, IAV


def zdmm() -> assembly.zdmm:
    return assembly.zdmm(resources.FASTA, resources.GFF3, HSV1(), IAV(), GRCm39())
