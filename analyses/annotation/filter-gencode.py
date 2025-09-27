import pickle

from biobit.core.loc import Interval
from intervaltree import IntervalTree
from joblib import Parallel, delayed
from pybedtools.cbedtools import defaultdict

import ld
import utils
from assemblies import CHM13v2, GRCm39


def job(assembly: str):
    obj = {"GRCm39": GRCm39, "CHM13v2": CHM13v2}[assembly]

    # Build an index of high confidence RefSeq lncRNAs
    high_confidence = defaultdict(IntervalTree)
    total_high_confidence, total_high_confidence_lncRNA = 0, 0
    for rna in obj.refseq.load().rnas.values():
        if rna.ind.startswith("rna-NR_"):
            total_high_confidence += 1
            if rna.attrs.biotype == "lnc_RNA":
                total_high_confidence_lncRNA += 1
            for exon in rna.exons:
                high_confidence[rna.loc.seqid, rna.loc.strand].addi(exon.start, exon.end)

    # Select only high confidence GENCODE lncRNAs
    allowed = set()
    autopassed = defaultdict(int)
    passed, discarded = 0, 0

    gencode = obj.gencode.load()
    for rna in gencode.rnas.values():
        if not ld.filters.is_well_defined(rna):
            continue

        if rna.attrs.type not in {
            "3prime_overlapping_ncRNA", "antisense", "bidirectional_promoter_lncRNA", "lincRNA", "macro_lncRNA",
            "non_coding", "processed_transcript", "sense_intronic", "sense_overlapping", "lncRNA",
        }:
            allowed.add(rna.ind)
            autopassed[rna.attrs.type] += 1
            continue

        # Check the overlap with high confidence RefSeq non-coding RNAs
        length, covered = sum(exon.len() for exon in rna.exons), 0
        for exon in rna.exons:
            maxcov = 0
            for overlap in high_confidence[rna.loc.seqid, rna.loc.strand].overlap(exon.start, exon.end):
                intersection = Interval(overlap.begin, overlap.end).intersection(exon)
                maxcov = max(maxcov, intersection.len() if intersection else 0)
            covered += maxcov

        if covered / length >= 0.75:
            allowed.add(rna.ind)
            passed += 1
        else:
            discarded += 1

    # Print a short report
    print(assembly)
    print(f"\tHigh confidence RefSeq non-coding RNAs: {total_high_confidence} (lncRNA: {total_high_confidence_lncRNA})")
    print(f"\tPassed: {passed}, Discarded: {discarded}")
    print("\tAutopassed:")
    for k, v in autopassed.items():
        print(f"\t\t{k}: {v}")

    saveto = ld.paths.gencode.pkl[assembly]
    saveto.parent.mkdir(parents=True, exist_ok=True)

    # Save as a pkl file
    gencode.rnas = {ind: rna for ind, rna in gencode.rnas.items() if ind in allowed}
    for gene in gencode.genes.values():
        transcripts = frozenset(ind for ind in gene.transcripts if ind in allowed)
        object.__setattr__(gene, "transcripts", transcripts)

    with open(saveto, "wb") as stream:
        pickle.dump(gencode, stream, protocol=pickle.HIGHEST_PROTOCOL)

    # Save as a BED file
    bed = utils.bed.blocks.from_iterable((
        (rna.loc.seqid, rna.exons, rna.loc.strand.to_orientation(), rna.attrs.name, None)
        for ind, rna in gencode.rnas.items()
    ))
    utils.bed.tbindex(bed, saveto.with_suffix(".bed.gz"))


Parallel(n_jobs=-1)(delayed(job)(assembly) for assembly in ["GRCm39", "CHM13v2"])
