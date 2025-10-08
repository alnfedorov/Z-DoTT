import pickle
from collections import defaultdict

from biobit.core.loc import Interval
from biobit.toolkit import annotome as at

from assemblies import CHM13v2, GRCm39


def hook(seqmap):
    def _hook(type: str, loc: at.transcriptome.Location, attributes: dict[str, str]):
        seqid = seqmap(loc.seqid)
        if not seqid:
            return None
        loc = at.transcriptome.Location(seqid, loc.strand, loc.start, loc.end)

        match type:
            case "CDS":
                attributes = {
                    "ID": attributes["ID"], "Parent": attributes["Parent"],
                    "product": attributes.get("product", None), "partial": attributes.get("partial", "false")
                }
            case "exon":
                # if "transcript_id" in attributes:
                #     assert attributes['Parent'] == f"rna-{attributes['transcript_id']}", attributes
                attributes = {"ID": attributes["ID"], "Parent": attributes["Parent"]}
            case "gene" | "pseudogene":
                assert attributes["gbkey"] == "Gene" and attributes["gene"] == attributes["Name"], attributes
                if attributes.pop("pseudo", None) == "true":
                    assert type == "pseudogene", attributes
                type = "gene"
            case "transcript":
                if attributes.pop("pseudo", None) == "true":
                    attributes["biotype"] = "pseudogene"
                else:
                    attributes["biotype"] = "unknown"
            case "miRNA":
                attributes["Parent"] = f"gene-{attributes['gene']}"
                attributes["biotype"] = "miRNA"
                type = "transcript"
            case "primary_transcript":
                attributes["biotype"] = "miRNA_primary_transcript"
                type = "transcript"
            case "mRNA" | "lnc_RNA" | "snRNA" | "snoRNA" | "antisense_RNA" | "tRNA" | "rRNA" | "telomerase_RNA" | \
                 "vault_RNA" | "Y_RNA" | "RNase_MRP_RNA" | "RNase_P_RNA" | "scRNA" | \
                 "V_gene_segment" | "C_gene_segment" | "J_gene_segment" | "D_gene_segment":
                if attributes.pop("pseudo", None) == "true":
                    attributes["biotype"] = f"{type}_pseudogene"
                else:
                    attributes["biotype"] = type
                type = "transcript"
            case "ncRNA":
                assert "Note" in attributes, attributes
                attributes["biotype"] = attributes.pop("Note")
                if "%3B" in attributes["biotype"]:
                    attributes["biotype"], _ = attributes["biotype"].split("%3B", 1)
                type = "transcript"
            case _:
                raise ValueError(f"Unknown type: {type}")

        if "Name" in attributes and "standard_name" in attributes:
            assert attributes["Name"] == attributes["standard_name"], attributes
            attributes.pop("standard_name")
        elif "standard_name" in attributes and "Name" not in attributes:
            attributes["Name"] = attributes.pop("standard_name")

        return type, loc, attributes

    return _hook


for assembly in GRCm39, CHM13v2:
    print(f"Processing {assembly.name} RefSeq GFF")
    records = at.preprocess_gff(
        assembly._refseq.gff,
        ignore_sources=set(),
        ignore_types={
            "biological_region", "enhancer", "silencer", "transcriptional_cis_regulatory_region",
            "protein_binding_site", "nucleotide_motif", "non_allelic_homologous_recombination_region",
            "recombination_feature", "promoter", "sequence_feature", "meiotic_recombination_region",
            "mobile_genetic_element", "DNaseI_hypersensitive_site", "conserved_region", "origin_of_replication",
            "tandem_repeat", "repeat_instability_region", "mitotic_recombination_region", "enhancer_blocking_element",
            "sequence_alteration", "TATA_box", "region", "response_element", "chromosome_breakpoint",
            "sequence_secondary_structure", "locus_control_region", "matrix_attachment_site",
            "epigenetically_modified_region", "replication_regulatory_region", "direct_repeat", "insulator",
            "minisatellite", "repeat_region", "CAAT_signal", "dispersed_repeat", "microsatellite", "inverted_repeat",
            "nucleotide_cleavage_site", "sequence_comparison", "GC_rich_promoter_region", "replication_start_site",
            "imprinting_control_region", "regulatory_region", "CAGE_cluster", "TSS", "sequence_alteration_artifact",
            "centromere", "match", "cDNA_match", "D_loop"
        },
        hook=hook(assembly.seqid.from_refseq), ind_key=lambda _, x: x["ID"]
    )
    for key in "CDS", "transcript", "gene":
        unique = set(x for matches in records[key].values() for _, x, _ in matches)
        print(f"\t{key} sources: {unique}")

    for name, key in ("transcript", "biotype"), ("gene", "gene_biotype"):
        unique = set(x[key] for matches in records[name].values() for _, _, x in matches)
        print(f"\t{name} biotypes: {unique}")

    # Parse CDS records
    cds, tid2cds = [], defaultdict(list)
    for ind, matches in records["CDS"].items():
        locations, sources, attributes = zip(*matches)

        assert len({(loc.seqid, loc.strand) for loc in locations}) == 1, locations
        seqid, strand = locations[0].seqid, locations[0].strand
        blocks = sorted(Interval(loc.start, loc.end) for loc in locations)

        assert len(set(sources)) == 1, sources
        source = sources[0]

        assert len(set(tuple(x.items()) for x in attributes)) == 1, attributes
        attributes = attributes[0]

        parents = set()
        for tid in attributes["Parent"].split(","):
            tid2cds[tid].append(ind)
            parents.add(tid)

        attrs = assembly._refseq.AttrCDS(source, bool(attributes["partial"]), attributes["product"], frozenset(parents))
        loc = at.transcriptome.Location(seqid, strand, blocks[0].start, blocks[-1].end)
        cds.append(at.transcriptome.CDS(ind, loc, attrs, tuple(blocks)))
    cds = at.transcriptome.CDSBundle(cds)

    # Parse exon records
    exons = defaultdict(list)
    for ind, matches in records["exon"].items():
        for loc, _, attributes in matches:
            for parent in attributes["Parent"].split(","):
                exons[parent].append(Interval(loc.start, loc.end))

    # Parse transcript records
    RNA, gid2tid = [], defaultdict(set)
    for ind, matches in records["transcript"].items():
        assert len(matches) == 1, matches
        loc, source, attributes = matches[0]

        rna_exons = sorted(exons.pop(ind), key=lambda x: x.start)
        rna_exons = tuple(x for x in rna_exons)

        parent = attributes.pop("Parent")
        assert "," not in parent, parent
        gid2tid[parent].add(ind)

        tags = [x for x in attributes.pop("tag", "").split(",") if x]
        attrs = assembly._refseq.AttrRNA(
            source, attributes.pop("Name", None), attributes.pop("product", None),
            bool(attributes.pop("partial", "false")), attributes.pop("biotype"), frozenset(tags),
            attributes.pop("experiment", None)
        )
        RNA.append(at.transcriptome.RNA(ind, loc, attrs, parent, rna_exons))
    RNA = at.transcriptome.RNABundle(RNA)

    # Parse gene records
    genes, singletons = [], set()
    for ind, matches in records["gene"].items():
        assert len(matches) == 1, matches
        location, source, attributes = matches[0]

        # Name, description, gene_biotype, gene_synonym, partial,
        synonyms = [x for x in attributes.pop("gene_synonym", "").split(",") if x]
        attrs = assembly._refseq.AttrGene(
            source, attributes.pop("Name"), attributes.pop("description", None), attributes.pop("gene_biotype"),
            bool(attributes.pop("partial", "false")), frozenset(synonyms)
        )
        transcripts = frozenset(gid2tid.pop(ind, []))
        if not transcripts:
            singletons.add(ind)
        genes.append(at.transcriptome.Gene(ind, location, attrs, transcripts))
    assert len(gid2tid) == 0, gid2tid
    genes = at.transcriptome.GeneBundle(genes)

    if singletons:
        print(f"Singleton genes (N={len(singletons)}): {singletons}")

    annotome = at.Annotome(assembly.name, "RefSeq", genes, RNA, cds)
    with open(assembly._refseq.index, "wb") as stream:
        pickle.dump(annotome, stream)
