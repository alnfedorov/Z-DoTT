from pathlib import Path
from typing import Literal

from attr import define, field
from biobit.toolkit import annotome as at

ROOT = Path(__file__).parent

type GeneSource = Literal['Liftoff', 'ENSEMBL', 'HAVANA']

type GeneType = Literal[
    'IG_C_gene', 'IG_C_pseudogene', 'IG_D_gene', 'IG_D_pseudogene', 'IG_J_gene', 'IG_J_pseudogene', 'IG_LV_gene',
    'IG_V_gene', 'IG_V_pseudogene', 'IG_pseudogene', 'Mt_rRNA', 'Mt_tRNA', 'TEC', 'TR_C_gene', 'TR_D_gene', 'TR_J_gene',
    'TR_J_pseudogene', 'TR_V_gene', 'TR_V_pseudogene', 'artifact', 'lncRNA', 'miRNA', 'misc_RNA',
    'processed_pseudogene', 'protein_coding', 'pseudogene', 'rRNA', 'rRNA_pseudogene', 'ribozyme', 'sRNA', 'scRNA',
    'scaRNA', 'snRNA', 'snoRNA', 'transcribed_processed_pseudogene', 'transcribed_unitary_pseudogene',
    'transcribed_unprocessed_pseudogene', 'translated_processed_pseudogene', 'translated_unprocessed_pseudogene',
    'unitary_pseudogene', 'unprocessed_pseudogene', 'vault_RNA'
]

type GeneLevel = Literal[1, 2, 3]


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrGene:
    # Liftoff attributes
    coverage: float
    sequence_identity: float
    copy_ind: int | None
    # GENCODE attributes
    source: GeneSource
    level: GeneLevel
    name: str
    type: GeneType

    def __attrs_post_init__(self):
        if self.coverage < 0 or self.coverage > 1:
            raise ValueError(f"Invalid coverage: {self.coverage}")
        if self.sequence_identity < 0 or self.sequence_identity > 1:
            raise ValueError(f"Invalid sequence identity: {self.sequence_identity}")
        if self.type not in GeneType.__value__.__args__:
            raise ValueError(f"Invalid biotype: {self.type}")
        if self.source not in GeneSource.__value__.__args__:
            raise ValueError(f"Invalid source: {self.source}")


type RNASource = Literal['Liftoff', 'ENSEMBL', 'HAVANA']

type RNAType = Literal[
    'IG_C_gene', 'IG_C_pseudogene', 'IG_D_gene', 'IG_D_pseudogene', 'IG_J_gene', 'IG_J_pseudogene', 'IG_LV_gene',
    'IG_V_gene', 'IG_V_pseudogene', 'IG_pseudogene', 'Mt_rRNA', 'Mt_tRNA', 'TEC', 'TR_C_gene', 'TR_D_gene', 'TR_J_gene',
    'TR_J_pseudogene', 'TR_V_gene', 'TR_V_pseudogene', 'artifact', 'lncRNA', 'miRNA', 'misc_RNA', 'non_stop_decay',
    'nonsense_mediated_decay', 'processed_pseudogene', 'processed_transcript', 'protein_coding',
    'protein_coding_CDS_not_defined', 'protein_coding_LoF', 'pseudogene', 'rRNA', 'rRNA_pseudogene', 'retained_intron',
    'ribozyme', 'sRNA', 'scRNA', 'scaRNA', 'snRNA', 'snoRNA', 'transcribed_processed_pseudogene',
    'transcribed_unitary_pseudogene', 'transcribed_unprocessed_pseudogene', 'translated_processed_pseudogene',
    'translated_unprocessed_pseudogene', 'unitary_pseudogene', 'unprocessed_pseudogene', 'vault_RNA'
]

type RNATag = Literal[
    '3_nested_supported_extension', '3_standard_supported_extension', '454_RNA_Seq_supported',
    '5_nested_supported_extension', '5_standard_supported_extension', 'CAGE_supported_TSS', 'CCDS', 'Ensembl canonical',
    'GENCODE basic', 'GENCODE primary', 'MANE Plus Clinical', 'MANE Select', 'NAGNAG_splice_site', 'NMD_exception',
    'NMD_likely_if_extended', 'RNA_Seq_supported_only', 'RNA_Seq_supported_partial', 'RP_supported_TIS', 'TAGENE',
    'alternative_3_UTR', 'alternative_5_UTR', 'appris_alternative_1', 'appris_alternative_2', 'appris_principal_1',
    'appris_principal_2', 'appris_principal_3', 'appris_principal_4', 'appris_principal_5', 'bicistronic', 'cds_end_NF',
    'cds_start_NF', 'dotter_confirmed', 'downstream_ATG', 'exp_conf', 'inferred_exon_combination',
    'inferred_transcript_model', 'low_sequence_quality', 'mRNA_end_NF', 'mRNA_start_NF', 'nested_454_RNA_Seq_supported',
    'non_ATG_start', 'non_canonical_TEC', 'non_canonical_U12', 'non_canonical_conserved',
    'non_canonical_genome_sequence_error', 'non_canonical_other', 'non_canonical_polymorphism',
    'non_submitted_evidence', 'not_best_in_genome_evidence', 'not_organism_supported', 'overlapping_uORF',
    'pseudo_consens', 'readthrough_transcript', 'retained_intron_CDS', 'retained_intron_final', 'retained_intron_first',
    'seleno', 'sequence_error', 'stop_codon_readthrough', 'upstream_ATG', 'upstream_uORF'
]

type RNATSL = Literal[1, 2, 3, 4, 5] | None
type RNALevel = Literal[1, 2, 3]


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrRNA:
    # Liftoff attributes
    copy_ind: int | None
    extra_copy_number: int
    # GENCODE attributes
    source: RNASource
    level: RNALevel
    name: str
    type: RNAType
    tags: frozenset[RNATag]
    TSL: RNATSL
    CDS: frozenset[str]

    def __attrs_post_init__(self):
        if self.source not in RNASource.__value__.__args__:
            raise ValueError(f"Invalid source: {self.source}")
        if self.type not in RNAType.__value__.__args__:
            raise ValueError(f"Invalid biotype: {self.type}")
        if any(x not in RNATag.__value__.__args__ for x in self.tags):
            raise ValueError(f"Invalid tags: {self.tags}")
        if self.TSL not in RNATSL.__value__.__args__:
            raise ValueError(f"Invalid TSL: {self.TSL}")


type CDSSource = Literal['Liftoff', 'ENSEMBL', 'HAVANA']


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrCDS:
    # Liftoff attributes
    copy_ind: int | None
    # GENCODE attributes
    source: CDSSource
    transcripts: frozenset[str] = field(converter=lambda x: frozenset(x))

    def __attrs_post_init__(self):
        if self.source not in CDSSource.__value__.__args__:
            raise ValueError(f"Invalid source: {self.source}")


type GencodeAnnotome = at.Annotome[AttrGene, AttrRNA, AttrCDS]
