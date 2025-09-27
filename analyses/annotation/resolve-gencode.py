import copy
import pickle
from collections import defaultdict
from typing import Any

import pybedtools
from biobit.core.loc import Interval, Orientation
from biobit.toolkit import countit
from joblib import Parallel, delayed

import ld
import utils
from assemblies import HSV1, IAV

REV_MAPPING = {
    'pseudogene': [
        'IG_D_pseudogene', 'translated_unprocessed_pseudogene', 'IG_C_pseudogene', 'TR_J_pseudogene',
        'IG_V_pseudogene', 'TR_V_pseudogene', 'transcribed_unitary_pseudogene',
        'transcribed_unprocessed_pseudogene', 'translated_processed_pseudogene', 'transcribed_processed_pseudogene',
        'unitary_pseudogene', 'unprocessed_pseudogene', 'processed_pseudogene', 'IG_pseudogene', 'IG_J_pseudogene',
        'rRNA_pseudogene'
    ],
    'mRNA': [
        'IG_D_gene', 'TR_D_gene', 'IG_J_gene', 'IG_V_gene', 'IG_C_gene', 'TR_V_gene', 'TR_J_gene', 'TR_C_gene',
        'IG_LV_gene', 'protein_coding',
    ],
    'mRNA defective': [
        'non_stop_decay', 'protein_coding_LoF', 'nonsense_mediated_decay',
        'retained_intron', 'protein_coding_CDS_not_defined',
    ],
    'smRNA': [
        'scRNA', 'vault_RNA', 'snoRNA', 'snRNA', 'miRNA', 'scaRNA', "sRNA"
    ],
    'lncRNA': [
        'lncRNA', 'processed_transcript'
    ]
}
MAPPING = {v: k for k, vs in REV_MAPPING.items() for v in vs}


def resolve_inner(keys: set[str]):
    initial = copy.deepcopy(keys)
    assert len(keys) > 0

    # Single biotype - no need to resolve
    if len(keys) == 1:
        return keys

    keys = {MAPPING.get(k, k) for k in keys if k != "intergenic"}

    # Singletons
    for k in "MT", "7SK", "7SL":
        if k in keys:
            return {k}

    # NMS for biotypes
    for core, suppressible in [
        ("mRNA", ["mRNA defective", "pseudogene"]),
        ("CDS", ["3'UTR", "5'UTR", "intron"]),
        ("mRNA defective", ["3'UTR", "5'UTR", "CDS", "pseudogene", "TEC", "intron"]),
        ("lncRNA", ["TEC", "mRNA defective", "pseudogene", "intron"]),
        ("pseudogene", ["mRNA defective", "TEC", "intron"]),
        ("ribozyme", ["lncRNA", "TEC", "pseudogene", "intron"]),
        ("TEC", ["intron"])
    ]:
        if core in keys:
            for k in suppressible:
                if k in keys:
                    keys.remove(k)

    # Resolve mRNAs
    if "mRNA" in keys:
        if "CDS" in keys:
            return {"CDS"}
        elif "5'UTR" in keys:
            return {"5'UTR"}
        elif "3'UTR" in keys:
            return {"3'UTR"}
        else:
            print(f"Unresolved mRNA: {initial} -> {keys}")
            return {"mRNA"}

    # Signletons
    for k in "smRNA", "rRNA", "tRNA", "misc_RNA":
        if k in keys:
            return {k}

    if len(keys) > 1:
        print(f"Unresolved: {initial} -> {keys}")
    return keys


def resolve(contig: str, orientation: Orientation, start: int, end: int, keys: set[str]) -> list[Any]:
    keys = resolve_inner(keys)
    return [(k, contig, orientation, start, end) for k in keys]


def job(assembly: str):
    assembly = utils.assembly.get(name=assembly)
    at = ld.load.gencode(assembly.name)

    # Construct the annotations
    annotation = defaultdict(lambda: defaultdict(list))

    # Gencode annotation
    for transcript in at.rnas.values():
        seqid, strand = transcript.loc.seqid, transcript.loc.strand.symbol()
        if transcript.attrs.name is not None and transcript.attrs.name.startswith("RN7SL"):
            annotation["7SL"][seqid, strand].append(transcript.exons[0])
        elif transcript.attrs.name is not None and transcript.attrs.name.startswith("RN7SK"):
            annotation["7SK"][seqid, strand].append(transcript.exons[0])
        else:
            # Add exons
            for ex in transcript.exons:
                annotation[transcript.attrs.type][seqid, strand].append(ex)

            # Add introns if available
            for prv, nxt in zip(transcript.exons[:-1], transcript.exons[1:]):
                if prv.end < nxt.start:
                    annotation["intron"][seqid, strand].append(Interval(prv.end, nxt.start))

            if transcript.attrs.type == "protein_coding" and not transcript.attrs.CDS:
                print(f"Missing CDS in the mRNA: {transcript}")

            # Add CDS & UTRs if available
            for CDS in transcript.attrs.CDS:
                cds = at.cds[CDS]
                assert cds.loc.seqid == seqid and cds.loc.strand.symbol() == strand
                for block in cds.blocks:
                    annotation["CDS"][seqid, strand].append(block)

                utr5 = []
                first, second = ("5'UTR", "3'UTR") if strand == "+" else ("3'UTR", "5'UTR")
                for ex in transcript.exons:
                    if ex.end < cds.loc.start:
                        utr5.append(ex)
                    else:
                        if ex.start < cds.loc.start < ex.end:
                            utr5.append(Interval(ex.start, cds.loc.start))
                        break
                annotation[first][seqid, strand].extend(utr5)

                utr3 = []
                for ex in transcript.exons[::-1]:
                    if ex.start > cds.loc.end:
                        utr3.append(ex)
                    else:
                        if ex.start < cds.loc.end < ex.end:
                            utr3.append(Interval(cds.loc.end, ex.end))
                        break
                annotation[second][seqid, strand].extend(utr3)

    seqsizes = assembly.seqid.sizes()

    # Mitochondria
    mt = seqsizes["chrM"]
    annotation["MT"]["chrM", "+"].append(Interval(0, mt))
    annotation["MT"]["chrM", "-"].append(Interval(0, mt))

    # Intergenic regions
    for contig, size in seqsizes.items():
        annotation["intergenic"][contig, "+"].append(Interval(0, size))
        annotation["intergenic"][contig, "-"].append(Interval(0, size))

    # Selected repmasker annotation
    for region in pybedtools.BedTool(assembly.repmasker):
        if region.name == "7SK":
            annotation["7SK"][region.chrom, region.strand].append(Interval(region.start, region.end))
        elif region.name == "7SLRNA":
            annotation["7SL"][region.chrom, region.strand].append(Interval(region.start, region.end))

    # Viral segments
    for virus in IAV, HSV1:
        for contig, size in virus.segments.items():
            for orient in Orientation.Forward, Orientation.Reverse:
                annotation[f"{virus.name}({orient})"][(contig, orient)].append(Interval(0, size))

    # Resolve the annotation
    annotation = {
        key: {
            loc: sorted(Interval.merge(segments)) for loc, segments in items.items()
        }
        for key, items in annotation.items()
    }
    annotation = countit.utils.resolve_annotation(annotation, resolve)

    # Save as a pickle
    saveto = ld.paths.resolved_annotation.pkl[assembly.name]
    saveto.parent.mkdir(parents=True, exist_ok=True)
    with open(saveto, "wb") as stream:
        pickle.dump(annotation, stream, protocol=pickle.HIGHEST_PROTOCOL)

    # Save as a bed file
    bed = pybedtools.BedTool([
        pybedtools.Interval(contig, s.start, s.end, strand=str(orient), name=name)
        for (name, *_), items in annotation.items()
        for (contig, orient), segments in items.items()
        for s in segments
    ])
    utils.bed.tbindex(bed, saveto.with_suffix(".bed.gz"))


Parallel(n_jobs=-1)(delayed(job)(assembly) for assembly in ("GRCm39", "CHM13v2"))
