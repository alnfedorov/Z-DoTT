from pathlib import Path

import pandas as pd
from BCBio import GFF
from Bio.SeqFeature import ExactPosition, FeatureLocation, SeqFeature

pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

ROOT = Path(__file__).parent

# Parse transcripts
file = ROOT / "supplements-transcripts.csv"
df = pd.read_csv(file)
df = df.rename(columns={"ID": "RNA-ID", "Num.Criteria": "Support"})

df['Note'] = [f"{name} [{js}]" if not pd.isna(js) else name
              for name, js in zip(df["Name"], df['Justification'])]
df['RNA-ID'] = df['RNA-ID'].apply(lambda x: x.replace("_RNA", "").replace("_", "-"))
# df['Gene-ID'] = df['Gene'].apply(lambda x: f"gene-{x}")
df['Gene-ID'] = df['RNA-ID'].apply(lambda x: f"gene-{x}")

df = df[['Gene-ID', 'RNA-ID', "Support", "Note", 'Location']]
assert len(df['RNA-ID'].unique()) == len(df['Gene-ID'].unique()) == len(df)

rnas, genes = set(df['RNA-ID']), set(df['Gene-ID'])
assert len(rnas & genes) == 0, rnas & genes

GENES = {k: [] for k in df['Gene-ID'].unique()}
TRANSCRIPTS = df

# Correct the GFF3 file
PROCESSED = set()
file = ROOT.joinpath("BK012101.1.gff3")
with open(file, 'r') as gff3:
    gff3 = list(GFF.parse(gff3))
    assert len(gff3) == 1
    gff3 = gff3[0]

OTHER_FEATURES = []
for gene in gff3.features:
    if gene.type in {"inverted_repeat", "sequence_feature", "stem_loop", "repeat_region", "region"}:
        OTHER_FEATURES.append(gene)
        continue
    assert gene.type == "gene", gene

    for RNA in gene.sub_features:
        if RNA.type == "CDS":
            continue

        assert RNA.type in {"mRNA", "ncRNA"}, RNA
        assert isinstance(RNA.location.start, ExactPosition) and isinstance(RNA.location.end, ExactPosition)

        start, end = RNA.location.start, RNA.location.end
        strand = "+" if RNA.location.strand > 0 else "-"

        exons = sorted([x for x in RNA.sub_features if x.type == 'exon'], key=lambda x: x.location.start)
        location = "|".join(f"{x.location.start}-{x.location.end}" for x in exons)
        location = f"JN555585{strand}:{location}"

        # Match with the transcript
        transcript = TRANSCRIPTS[TRANSCRIPTS['Location'] == location]
        if len(transcript) == 0:
            # Strange record without a matched transcript from the paper
            continue
        assert len(transcript) == 1, transcript

        transcript = transcript.iloc[0]
        RNA.qualifiers['Note'] = transcript['Note']
        RNA.qualifiers['support'] = transcript['Support']

        # - Fix IDs
        RNA.id = transcript['RNA-ID']
        RNA.qualifiers['ID'] = transcript['RNA-ID']
        RNA.qualifiers['Parent'] = transcript['Gene-ID']
        RNA.qualifiers['gene'] = transcript['Gene-ID']
        for x in RNA.sub_features:
            assert len(x.sub_features) == 0
            x.qualifiers['Parent'] = RNA.qualifiers['ID']
            x.qualifiers['gene'] = transcript['Gene-ID']

        # Drop ORFs
        RNA.sub_features = [x for x in RNA.sub_features if x.type != 'CDS']

        PROCESSED.add(transcript['RNA-ID'])
        GENES[transcript['Gene-ID']].append(RNA)

assert len(PROCESSED) == len(TRANSCRIPTS), TRANSCRIPTS[~TRANSCRIPTS['RNA-ID'].isin(PROCESSED)]

# Infer genes
genes = []
for geneid, transcript in GENES.items():
    assert all(x.location.strand == transcript[0].location.strand for x in transcript)
    strand = transcript[0].location.strand
    start, end = min(transcript, key=lambda x: x.location.start).location.start, \
        max(transcript, key=lambda x: x.location.end).location.end
    name = geneid.replace("gene-", '')
    gene = SeqFeature(
        FeatureLocation(start, end, strand), "gene",
        qualifiers={"ID": geneid, "Name": name, "gene": geneid, "source": ['tpg']}
    )
    gene.sub_features = transcript
    genes.append(gene)
gff3.features = OTHER_FEATURES + genes

saveto = ROOT.joinpath("annotation.gff3")
with open(saveto, "w") as out_handle:
    GFF.write([gff3], out_handle)
