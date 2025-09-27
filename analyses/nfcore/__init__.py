from pathlib import Path

# Internal datasets
from .rnaseq.B256178 import B256178  # IAV infection of MEFs (batch 1)
from .rnaseq.B261790 import B261790  # HSV-1 infection of MEFs (batch 1, Z22 RIP-seq)
from .rnaseq.B319096 import B319096  # HSV-1 infection of MEFs (batch 2, FLAG RIP-seq)
from .rnaseq.B831009 import B831009  # HSV-1/IAV RIP from HT-29s (FLAG + Z22 RIP-seq)

# SRA Projects
from .rnaseq.PRJEB75711 import PRJEB75711
from .rnaseq.PRJNA256013 import PRJNA256013
from .rnaseq.PRJNA382632 import PRJNA382632
from .rnaseq.PRJNA637636 import PRJNA637636

ROOT = Path(__file__).parent

BY_ASSEMBLY = {
    "CHM13v2": [
        PRJNA256013, PRJNA382632, PRJNA637636, PRJEB75711, B831009,
    ],
    "GRCm39": [
        B256178, B261790, B319096,
    ],
}
