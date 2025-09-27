import pickle

from joblib import Parallel, delayed
from pybedtools import BedTool

import ld
import utils.repeto
from stories import annotation


def derive_insulators(config: ld.Config):
    insulators = [(p.chrom, (p.start, p.end), p.strand) for p in BedTool(config.peaks.not_covered)]

    # Load all strong natural RNA "insulators" - splicing sites inside canonical mRNAs
    assembly = utils.assembly.get(organism=config.host)
    gencode = assembly.gencode.load()
    for rna in gencode.rnas.values():
        if annotation.filters.is_primary(rna) and rna.attrs.type == "protein_coding":
            for prv, nxt in zip(rna.exons[:-1], rna.exons[1:]):
                if prv.end < nxt.start:
                    insulators.append((rna.loc.seqid, (prv.end, prv.end + 1), rna.loc.strand.symbol()))
                    insulators.append((rna.loc.seqid, (nxt.start - 1, nxt.start), rna.loc.strand.symbol()))

    # Load curated insulators
    if config.dsRNA.insulators.exists():
        for p in BedTool(config.dsRNA.insulators):
            insulators.append((p.chrom, (p.start, p.end), "+"))
            insulators.append((p.chrom, (p.start, p.end), "-"))
    return config.ind, insulators


result = Parallel(n_jobs=-1)(delayed(derive_insulators)(config) for config in ld.Config.load())
result = dict(result)

ld.INSULATORS_CACHE.parent.mkdir(parents=True, exist_ok=True)
with open(ld.INSULATORS_CACHE, "wb") as f:
    pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)
