import pybedtools
from pybedtools import BedTool

import ld
from utils import bed

dsRNA = "Z-RNA[GRCm39]"

records = {
    "exclude": """
    """,
    "insulators": """
    """,
    "connectors": """
    """,
    "clusters": """
    """,
}

for name, allnew in records.items():
    fname = ld.RESULTS / dsRNA / "curated" / f"{name}.bed"
    coordinates = [bed.parse_coords(x) for x in allnew.split() if x]
    for c in coordinates:
        c.end += 1
    if fname.exists():
        coordinates += [pybedtools.Interval(x.chrom, x.start, x.end) for x in pybedtools.BedTool(fname)]
    coordinates = BedTool(coordinates).sort().merge().sort()

    with open(fname, 'w') as stream:
        for interval in coordinates:
            stream.write(f"{interval.chrom}\t{interval.start}\t{interval.end}\n")
