from collections import defaultdict
from pathlib import Path

from biobit.core.loc import Interval
from biobit.toolkit.annotome import Annotome
from biobit.toolkit.countit.utils import resolve_annotation
from intervaltree import IntervalTree
from pybedtools import BedTool


def annotate(assembly, anno: Annotome, sites: Path):
    # Fetch the genomic annotation
    regions = defaultdict(lambda: defaultdict(list))
    for rna in anno.rnas.values():
        # Genic
        for exon in rna.exons:
            regions["exon"][rna.loc.seqid, str(rna.loc.strand)].append(exon)
        for prv, nxt in zip(rna.exons[:-1], rna.exons[1:]):
            if prv.end < nxt.start:
                regions["intron"][rna.loc.seqid, str(rna.loc.strand)].append(Interval(prv.end, nxt.start))
        # Intergenic
        tes = rna.loc.tes
        match rna.loc.strand:
            case "+":
                regions["TES+1kb"][rna.loc.seqid, "+"].append(Interval(tes, tes + 1000))
                regions["TES+5kb"][rna.loc.seqid, "+"].append(Interval(tes, tes + 5000))
                regions["TES+10kb"][rna.loc.seqid, "+"].append(Interval(tes, tes + 10000))
            case "-":
                regions["TES+1kb"][rna.loc.seqid, "-"].append(Interval(tes - 1000, tes))
                regions["TES+5kb"][rna.loc.seqid, "-"].append(Interval(tes - 5000, tes))
                regions["TES+10kb"][rna.loc.seqid, "-"].append(Interval(tes - 10000, tes))
            case _:
                raise ValueError(f"Invalid strand: {rna.loc.strand}")

    # Mitochondria
    mt = assembly.seqid.sizes()["chrM"]
    regions["MT"]["chrM", "+"].append(Interval(0, mt))
    regions["MT"]["chrM", "-"].append(Interval(0, mt))

    # Resolve the annotation
    def resolution(_contig, _orientation, _start, _end, keys):
        for key in "MT", "exon", "intron", "TES+1kb", "TES+5kb", "TES+10kb":
            if key in keys:
                return [key]
        return ["intergenic"]

    resolved = resolve_annotation(regions, resolution)

    trees = defaultdict(IntervalTree)
    for key, anno in resolved.items():
        for (contig, orientation), rregions in anno.items():
            tkey = (contig, orientation.symbol())
            for region in rregions:
                trees[tkey].addi(region.start, region.end, data=key)

    # Annotate sites
    categories, mapping = defaultdict(list), {}
    for site in BedTool(sites):
        overlap = trees[site.chrom, site.strand].overlap(site.start, site.end)
        assert len(overlap) <= 1, (overlap, site.chrom, site.strand, site.start, site.end)

        cat = overlap.pop().data if overlap else "intergenic"
        categories[cat].append(site)
        mapping[site.chrom, site.strand, site.start] = cat

    return categories, mapping
