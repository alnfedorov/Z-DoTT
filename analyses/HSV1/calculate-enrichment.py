import pickle
from itertools import chain

from biobit.core import ngs
from biobit.toolkit import nfcore
from joblib import Parallel, delayed, cpu_count

import ld
from assemblies import HSV1
from stories import normalization
from stories.RIP.pcalling.ld import Config

ld.ENRICHMENT.mkdir(parents=True, exist_ok=True)

SCALING = normalization.median_of_ratios()
COMPARISONS = [x for x in Config.load() if 'Herpes simplex virus 1' in x.organism]


def job(cmp: Config):
    saveto = ld.ENRICHMENT / f"{cmp.ind}.pkl"
    if saveto.exists():
        return

    # 1. Extract path to BAM files
    _signal = [nfcore.rnaseq.extract.bam(x) for x in cmp.signal]
    _control = [nfcore.rnaseq.extract.bam(x) for x in cmp.control]

    assert all(
        isinstance(x[1], ngs.Layout.Paired) and
        x[1].strandedness == ngs.Strandedness.Reverse and x[1].orientation == ngs.MatesOrientation.Inward
        for x in chain(_signal, _control)
    )
    signal = [x[0].filename for x in _signal]
    control = [x[0].filename for x in _control]

    # 2. Calculate the scaling factor
    _, _, scfactors = cmp.scaling(SCALING[cmp.assembly])
    scfactors = {k: 1 / v for k, v in scfactors.items()}
    print(cmp.ind, scfactors)

    # 3. Calculate the enrichment
    enrichment = ld.enrichment.calculate(
        HSV1.segments, cmp.ind, ld.ENRICHMENT,
        signal, scfactors["signal"],
        control, scfactors["control"],
        cpus=cpu_count(), exclude_flags=2572, include_flags=3, minmapq=0
    )

    # 4. Save the enrichment
    with open(saveto, "wb") as stream:
        pickle.dump(enrichment, stream)


Parallel(n_jobs=-1)(delayed(job)(cmp) for cmp in COMPARISONS)
