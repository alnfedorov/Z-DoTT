import pickle
from collections import defaultdict

import ld
import utils
from stories import annotation

RNA = {}
for assembly in ld.comparisons.series:
    RNA[assembly] = defaultdict(list)
    gencode = utils.assembly.get(name=assembly).gencode.load()

    for (seqid, strand), cores in annotation.load.rna_cores(assembly).items():
        allowed = []
        for rna in cores:
            if gencode.genes[rna.gid].attrs.type != "protein_coding" or \
                    sum(exon.len() for exon in rna.exons) < ld.thr.min_ref_length:
                continue
            allowed.append(rna)
        RNA[assembly][seqid].extend(allowed)

RNA = {k: dict(v) for k, v in RNA.items()}

ld.RNA.parent.mkdir(parents=True, exist_ok=True)
with open(ld.RNA, 'wb') as stream:
    pickle.dump(RNA, stream)
