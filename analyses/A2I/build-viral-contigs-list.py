from pybedtools import Interval, BedTool

import ld
from assemblies import HSV1, IAV
from utils import fasta

intervals = []
for virus in HSV1, IAV:
    for contig, length in fasta.contigs(virus.fasta).items():
        intervals.append(Interval(contig, 0, length, strand="+"))
        intervals.append(Interval(contig, 0, length, strand="-"))

BedTool(intervals).sort().saveas(ld.c.VIRAL)
