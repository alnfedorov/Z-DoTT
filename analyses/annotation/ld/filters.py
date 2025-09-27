from biobit.toolkit.annotome.transcriptome import RNA

from assemblies import CHM13v2, GRCm39


def is_well_defined(rna: RNA) -> bool:
    assert isinstance(rna.attrs, (CHM13v2.gencode.AttrRNA, GRCm39.gencode.AttrRNA))
    allow = len(rna.attrs.tags & {
        "retained_intron_CDS", "retained_intron_final", "retained_intron_first", "readthrough_transcript",
    }) == 0
    allow &= rna.attrs.type not in {
        'retained_intron', 'TEC',
        # Pseudogenes (even though some of them might be expressed - they are undetectable in our data)
        'IG_C_pseudogene', 'IG_D_pseudogene', 'IG_J_pseudogene', 'IG_V_pseudogene',
        'IG_pseudogene', 'TR_J_pseudogene', 'TR_V_pseudogene', 'processed_pseudogene', 'pseudogene', 'rRNA_pseudogene',
        'transcribed_processed_pseudogene', 'transcribed_unitary_pseudogene', 'transcribed_unprocessed_pseudogene',
        'translated_processed_pseudogene', 'translated_unprocessed_pseudogene', 'unitary_pseudogene',
        'unprocessed_pseudogene'
    }
    # Drop automatic annotation
    allow &= rna.attrs.level <= 2
    # Keep only GENCODE basic transcripts
    allow &= "GENCODE basic" in rna.attrs.tags

    # Allow manually verified transcripts
    verified = rna.ind in {
        # U1 snRNAs
        "ENSMUST00000239775.1", "ENSMUST00001239542.1",
    }
    return allow or verified


def is_primary(rna: RNA) -> bool:
    assert isinstance(rna.attrs, (CHM13v2.gencode.AttrRNA, GRCm39.gencode.AttrRNA))
    return is_canonical(rna) or len(rna.attrs.tags & {"GENCODE primary", "MANE Plus Clinical"}) != 0


def is_canonical(rna: RNA) -> bool:
    assert isinstance(rna.attrs, (CHM13v2.gencode.AttrRNA, GRCm39.gencode.AttrRNA))
    return len(rna.attrs.tags & {"Ensembl canonical", "MANE Select"}) != 0


def is_transcription_boundary(rna: RNA) -> bool:
    # Small RNAs are not captured efficiently by the 150bp PE RNA-seq
    if rna.loc.end - rna.loc.start < 100:
        return False

    if isinstance(rna.attrs, (GRCm39.refseq.AttrRNA, CHM13v2.refseq.AttrRNA)):
        if rna.attrs.biotype in {
            # Artefacts
            'unknown',
            # Pseudogenes
            'pseudogene', 'lnc_RNA_pseudogene', 'V_gene_segment_pseudogene', 'C_gene_segment_pseudogene',
            'J_gene_segment_pseudogene', 'D_gene_segment_pseudogene',
            # rRNAs are depleted prior to sequencing and during mapping
            'rRNA',
        }:
            return False
    elif isinstance(rna.attrs, (CHM13v2.gencode.AttrRNA, GRCm39.gencode.AttrRNA)):
        if not is_well_defined(rna) or rna.attrs.type in {
            # rRNA
            'rRNA', 'Mt_rRNA',
        }:
            return False
    return True
