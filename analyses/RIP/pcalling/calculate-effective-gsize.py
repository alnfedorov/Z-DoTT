import gzip
import pickle

import khmer
from Bio import SeqIO
from joblib import Parallel, delayed

import ld
from assemblies import GRCm39, CHM13v2

ERROR_RATE = 0.01
KSIZE = 150


def job(assembly: str, contig: SeqIO.SeqRecord):
    counter = khmer.HLLCounter(ERROR_RATE, KSIZE)
    counter.consume_string(str(contig.seq))

    cardinality = counter.estimate_cardinality()
    length = len(contig.seq)
    notn = length - contig.seq.count("N") - contig.seq.count("n")

    print(f"{contig.id}")
    print(f"\tUnique {KSIZE}-mers: {cardinality} ({cardinality / length * 100:.2f}%)")
    print(f"\tNot-N: {notn} ({notn / length * 100:.2f}%)")
    es = min(cardinality, notn)
    assert es <= length
    print(f"\tFinal effective size = {es} ({es / length * 100:.2f}%)")
    return assembly, contig.id, es


workload = []
for assembly in GRCm39, CHM13v2:
    with gzip.open(assembly.fasta, 'rt') as stream:
        for contig in SeqIO.parse(stream, 'fasta'):
            workload.append((assembly.fasta, contig))

result = Parallel(n_jobs=-1)(
    delayed(job)(assembly, contig) for assembly, contig in workload
)

effective_size = {GRCm39.name: {}, CHM13v2.name: {}}
for assembly, contig, es in result:
    effective_size[assembly][contig] = es

ld.EFFECTIVE_GENOME_SIZE.parent.mkdir(parents=True, exist_ok=True)
with open(ld.EFFECTIVE_GENOME_SIZE, 'wb') as stream:
    pickle.dump(effective_size, stream)
