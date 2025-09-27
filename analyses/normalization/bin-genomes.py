import pickle
from collections import defaultdict

from biobit.core.loc import Interval, Strand
from pybedtools import BedTool

import ld
import utils
from assemblies import IAV, CHM13v2, GRCm39, HSV1
from stories import annotation
from stories.normalization.ld import NormBin

ld.RESULTS.mkdir(parents=True, exist_ok=True)


def bin_genome(assembly) -> dict[tuple[str, Strand], tuple[NormBin, ...]]:
    gencode = annotation.load.gencode(assembly.name)

    # Construct the list of introns/exons for well-defined genes
    genes = defaultdict(list)
    for rna in gencode.rnas.values():
        if rna.loc.seqid == 'chrM':
            continue
        for exon in rna.exons:
            genes["exon", rna.gene].append(exon)
        for prv, nxt in zip(rna.exons[:-1], rna.exons[1:]):
            if nxt.start != prv.end:
                genes["intron", rna.gene].append(Interval(prv.end, nxt.start))

    # Merge the annotated regions
    bins = defaultdict(list)
    for (cat, gene), intervals in genes.items():
        record = gencode.genes[gene]
        seqid, strand = record.loc.seqid, record.loc.strand

        intervals = Interval.merge(intervals)
        bins[seqid, strand].append(NormBin(f"{gene}[{cat}]", cat, seqid, strand, tuple(intervals)))

    # Derive intergenic regions
    seqsizes = assembly.seqid.sizes()
    for seqid, seqlen in seqsizes.items():
        if seqid == 'chrM':
            continue

        for strand in Strand.Forward, Strand.Reverse:
            annotated = Interval.merge([
                interval for partition in bins.get((seqid, strand), []) for interval in partition.intervals
            ])

            if len(annotated) == 0:
                bins[seqid, strand].append(NormBin(
                    f"{seqid}:{0}-{seqlen}", "intergenic", seqid, strand, (Interval(0, seqlen),)
                ))
                continue

            # Intergenic region before the first interval
            if annotated[0].start != 0:
                bins[seqid, strand].append(NormBin(
                    f"{seqid}:{0}-{annotated[0].start}", "intergenic", seqid, strand, (Interval(0, annotated[0].start),)
                ))

            # Intergenic regions between annotated intervals
            for prv, nxt in zip(annotated[:-1], annotated[1:]):
                assert prv.end < nxt.start
                if nxt.start != prv.end:
                    bins[seqid, strand].append(NormBin(
                        f"{seqid}:{prv.end}-{nxt.start}", "intergenic", seqid, strand, (Interval(prv.end, nxt.start),)
                    ))

            # Intergenic region post the last interval
            if annotated[-1].end != seqlen:
                bins[seqid, strand].append(NormBin(
                    f"{seqid}:{annotated[-1].end}-{seqlen}", "intergenic", seqid, strand,
                    (Interval(annotated[-1].end, seqlen),)
                ))

    for strand in Strand.Forward, Strand.Reverse:
        # Viruses
        for virus in IAV, HSV1:
            for seqid, seqlen in virus.segments.items():
                bins[seqid, strand].append(NormBin(
                    f"{virus.name}[{strand}]", 'vRNA', seqid, strand, (Interval(0, seqlen),)
                ))

        # Mitochondrial genome
        bins["chrM", strand].append(NormBin(
            f"chrM[{strand}]", "MT", "chrM", strand, (Interval(0, seqsizes["chrM"]),)
        ))

    # Save generated bins
    pybed = []
    for parts in bins.values():
        for element in parts:
            pybed.append(utils.bed.blocks.make(
                element.seqid, element.intervals, element.strand.symbol(), element.ind
            ))
    utils.bed.tbindex(BedTool(pybed).sort(), ld.RESULTS / f"{assembly.name}-bin.bed.gz")

    return {k: tuple(v) for k, v in bins.items()}


bins = {}
for assembly in ld.SERIES:
    assembly = {"GRCm39": GRCm39, "CHM13v2": CHM13v2}[assembly]
    bins[assembly.name] = bin_genome(assembly)

ld.BINS.parent.mkdir(parents=True, exist_ok=True)
with open(ld.BINS, 'wb') as stream:
    pickle.dump(bins, stream)
