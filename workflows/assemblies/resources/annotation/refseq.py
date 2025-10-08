from pathlib import Path
from typing import Literal

from attr import define, field
from biobit.toolkit import annotome as at

ROOT = Path(__file__).parent

type GeneSource = Literal[
    'BestRefSeq', 'BestRefSeq%2CGnomon', 'BestRefSeq%2Ccmsearch', 'Curated Genomic', 'Curated Genomic%2Ccmsearch',
    'Gnomon', 'RefSeq', 'cmsearch', 'tRNAscan-SE'
]

type GeneBiotype = Literal[
    'C_region', 'C_region_pseudogene', 'D_segment', 'D_segment_pseudogene', 'J_segment', 'J_segment_pseudogene',
    'RNase_MRP_RNA', 'RNase_P_RNA', 'V_segment', 'V_segment_pseudogene', 'Y_RNA', 'antisense_RNA', 'lncRNA', 'miRNA',
    'misc_RNA', 'ncRNA', 'ncRNA_pseudogene', 'other', 'protein_coding', 'pseudogene', 'rRNA', 'scRNA', 'snRNA',
    'snoRNA', 'tRNA', 'telomerase_RNA', 'transcribed_pseudogene', 'vault_RNA'
]


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrGene:
    source: GeneSource
    name: str
    description: str | None
    biotype: GeneBiotype
    partial: bool
    synonyms: frozenset[str] = field(converter=lambda x: frozenset(x))

    def __attrs_post_init__(self):
        if self.source not in GeneSource.__value__.__args__:
            raise ValueError(f"Invalid source: {self.source}")
        if self.biotype not in GeneBiotype.__value__.__args__:
            raise ValueError(f"Invalid biotype: {self.biotype}")


type RNASource = Literal['BestRefSeq', 'Curated Genomic', 'Gnomon', 'RefSeq', 'cmsearch', 'tRNAscan-SE']

type RNABiotype = Literal[
    'C_gene_segment', 'C_gene_segment_pseudogene', 'D_gene_segment', 'D_gene_segment_pseudogene', 'J_gene_segment',
    'J_gene_segment_pseudogene', 'RNase_MRP_RNA', 'RNase_P_RNA', 'V_gene_segment', 'V_gene_segment_pseudogene',
    'Y_RNA', 'antisense_RNA', 'lnc_RNA', 'lnc_RNA_pseudogene', 'mRNA', 'miRNA', 'miRNA_primary_transcript',
    'pseudogene', 'rRNA', 'scRNA', 'scaRNA', 'snRNA', 'snoRNA', 'tRNA', 'telomerase_RNA', 'unknown', 'vault_RNA'
]

type RNAExperiment = Literal[
    'COORDINATES: cap analysis [ECO:0007248]',
    'COORDINATES: cap analysis [ECO:0007248] and polyA evidence [ECO:0006239]',
    'COORDINATES: polyA evidence [ECO:0006239]'
]

type RNATag = Literal['MANE Plus Clinical', 'MANE Select', 'RefSeq Plus Clinical', 'RefSeq Select']


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrRNA:
    source: RNASource
    name: str | None
    product: str | None
    partial: bool
    biotype: RNABiotype
    tags: frozenset[RNATag]
    experiment: RNAExperiment | None

    def __attrs_post_init__(self):
        if self.source not in RNASource.__value__.__args__:
            raise ValueError(f"Invalid source: {self.source}")
        if self.biotype not in RNABiotype.__value__.__args__:
            raise ValueError(f"Invalid biotype: {self.biotype}")
        if any(x not in RNATag.__value__.__args__ for x in self.tags):
            raise ValueError(f"Invalid tags: {self.tags}")
        if self.experiment and self.experiment not in RNAExperiment.__value__.__args__:
            raise ValueError(f"Invalid experiment: {self.experiment}")


type CDSSource = Literal['BestRefSeq', 'Curated Genomic', 'Gnomon', 'RefSeq']


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrCDS:
    source: CDSSource
    partial: bool
    product: str | None
    transcripts: frozenset[str] = field(converter=lambda x: frozenset(x))

    def __attrs_post_init__(self):
        if self.source not in CDSSource.__value__.__args__:
            raise ValueError(f"Invalid source: {self.source}")


type RefSeqAnnotome = at.Annotome[AttrGene, AttrRNA, AttrCDS]
