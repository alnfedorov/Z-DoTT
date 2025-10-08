import pickle
from collections import defaultdict
from typing import Iterable

from biobit.core.loc import Interval
from biobit.toolkit import annotome as at

from assemblies import CHM13v2


def assert_attributes_match(attributes: Iterable[dict[str, str]]):
    unique = defaultdict(lambda: defaultdict(int))
    elements = 0

    failed = False
    for attrs in attributes:
        elements += 1
        for key, value in attrs.items():
            unique[key][value] += 1
    for key, values in unique.items():
        if len(values) > 1:
            print(f"WARNING: Inconsistent value for {key}: {values}")
            failed = True
        for value, count in values.items():
            if count != elements:
                print(f"WARNING: Inconsistent count for {key}={value}: {count}/{elements}")
                failed = True
    assert not failed, "Inconsistent attributes"


def hook(type: str, loc: at.transcriptome.Location, attributes: dict[str, str]):
    match type:
        case "CDS":
            for key in "exon_number", "exon_id":
                attributes.pop(key)
        case "exon":
            pass
        case "gene" | "ncRNA_gene" | "pseudogene":
            type = "gene"
        case "mRNA" | "lnc_RNA" | "transcript" | "pseudogenic_transcript" | \
             "unconfirmed_transcript" | "processed_transcript" | "gene_segment" | \
             "V_gene_segment" | "J_gene_segment" | "C_gene_segment" | "D_gene_segment" | \
             "miRNA" | "snoRNA" | "snRNA" | "rRNA" | "scRNA" | "tRNA":
            type = "transcript"
        case _:
            raise ValueError(f"Unknown type: {type}")
    return type, loc, attributes


assembly = CHM13v2
print(f"Processing {assembly.name} GENCODE GFF3")
records = at.preprocess_gff(
    assembly._gencode.gff3,
    ignore_sources={"GRCm39", "GRCh38", "cpg", "Eponine"},
    ignore_types={
        "five_prime_UTR", "three_prime_UTR", "biological_region", "chromosome", "scaffold", "start_codon", "stop_codon",
        "stop_codon_redefined_as_selenocysteine"
    },
    hook=hook, ind_key=lambda _, x: x["ID"]
)

sources = set(x for matches in records["CDS"].values() for _, x, _ in matches)
print(f"\tCDS sources: {sources}")

# Parse CDS records
cds, tid2cds = [], defaultdict(list)
for ind, matches in records["CDS"].items():
    locations, sources, attributes = zip(*matches)

    assert len({(loc.seqid, loc.strand) for loc in locations}) == 1, locations
    seqid, strand = locations[0].seqid, locations[0].strand
    blocks = sorted(Interval(loc.start, loc.end) for loc in locations)

    assert len(set(sources)) == 1, sources
    source = sources[0]

    assert_attributes_match(attributes)
    attributes = attributes[0]

    parents = set()
    for tid in attributes["Parent"].split(","):
        assert not tid.startswith("transcript:"), tid

        tid2cds[tid].append(ind)
        parents.add(tid)

    if "_" not in ind:
        copyind = None
    else:
        _, copyind = ind.split("_", maxsplit=1)
        copyind = int(copyind)

    attrs = assembly._gencode.AttrCDS(copyind, source, frozenset(parents))
    loc = at.transcriptome.Location(seqid, strand, blocks[0].start, blocks[-1].end)
    cds.append(at.transcriptome.CDS(ind, loc, attrs, tuple(blocks)))
cds = at.transcriptome.CDSBundle(cds)

# Parse exon records
exons = defaultdict(list)
for ind, matches in records["exon"].items():
    for location, _, attributes in matches:
        rank = int(attributes.pop("exon_number"))
        for parent in attributes["Parent"].split(","):
            assert not parent.startswith("transcript:"), parent
            exons[parent].append((rank, Interval(location.start, location.end)))

# Parse transcript records
ttypes = {x[0][2]['transcript_type'] for x in records["transcript"].values()}
print(f"\tTranscript types: {ttypes}")

tags = {tname for x in records["transcript"].values() for tname in x[0][2].get('tag', '').split(",") if tname}
print(f"\tTranscript tags: {tags}")

sources = set(x for matches in records["transcript"].values() for _, x, _ in matches)
print(f"\tTranscript sources: {sources}")

RNA, gid2tid = [], defaultdict(set)
for ind, matches in records["transcript"].items():
    assert len(matches) == 1, matches
    location, source, attributes = matches[0]

    rna_exons = sorted(exons.pop(ind), key=lambda x: x[0], reverse=(location.strand == "-"))
    for i in range(1, len(rna_exons)):
        # if rna_exons[i - 1][1].end >= rna_exons[i][1].start:
        #     print(f"WARNING: Overlapping exons in {ind}: {rna_exons[i - 1]} and {rna_exons[i]}")
        assert rna_exons[i - 1][1].end <= rna_exons[i][1].start, (rna_exons[i - 1], rna_exons[i])
    rna_exons = tuple(x[1] for x in rna_exons)

    parent = attributes.pop("Parent")
    assert not parent.startswith("gene:") and "," not in parent, parent
    gid2tid[parent].add(ind)

    # Parse the TSL
    tsl = attributes.pop("transcript_support_level", "NA")
    if " " in tsl:
        tsl, postfix = tsl.split(" ", maxsplit=1)
        assert postfix.startswith("(assigned to previous version") or len(postfix) == 0, (tsl, postfix)
    tsl = {"NA": None, "1": 1, "2": 2, "3": 3, "4": 4, "5": 5}[tsl]

    # Parse tags
    tags = [{"basic": "GENCODE basic", "Ensembl_canonical": "Ensembl canonical", "MANE_Select": "MANE Select",
             "GENCODE_Primary": "GENCODE primary", "MANE_Plus_Clinical": "MANE Plus Clinical"
             }.get(t, t) for t in attributes.pop("tag", "").split(",") if t]

    if "_" not in ind:
        copyind = None
    else:
        _, copyind = ind.split("_", maxsplit=1)
        copyind = int(copyind)

    attrs = assembly._gencode.AttrRNA(
        copyind, int(attributes.pop("extra_copy_number")), source, int(attributes.pop("level")),
        attributes.pop("transcript_name"), attributes.pop("transcript_type"), frozenset(tags), tsl,
        frozenset(tid2cds.pop(ind, ()))
    )
    RNA.append(at.transcriptome.RNA(ind, location, attrs, parent, tuple(rna_exons)))
assert len(exons) == 0, exons
assert len(tid2cds) == 0, tid2cds
RNA = at.annotome.RNABundle(RNA)

# Parse gene records
gtype = {x[0][2]['gene_type'] for x in records["gene"].values()}
print(f"\tGene biotypes: {gtype}")

sources = set(x for matches in records["gene"].values() for _, x, _ in matches)
print(f"\tGene sources: {sources}")

genes = []
for ind, matches in records["gene"].items():
    assert len(matches) == 1, matches
    location, source, attributes = matches[0]

    if "_" not in ind:
        copyind = None
    else:
        _, copyind = ind.split("_", maxsplit=1)
        copyind = int(copyind)

    attrs = assembly._gencode.AttrGene(
        float(attributes.pop("coverage")), float(attributes.pop("sequence_ID")),
        copyind, source, int(attributes.pop("level")), attributes.pop("gene_name"), attributes.pop("gene_type")
    )
    genes.append(at.transcriptome.Gene(ind, location, attrs, frozenset(gid2tid.pop(ind))))
assert len(gid2tid) == 0, gid2tid
genes = at.transcriptome.GeneBundle(genes)

annotome = at.Annotome(assembly.name, "GENCODE v47", genes, RNA, cds)
with open(assembly._gencode.index, "wb") as stream:
    pickle.dump(annotome, stream)
