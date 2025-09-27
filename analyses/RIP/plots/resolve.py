from itertools import chain

import pandas as pd

import ld
from assemblies import CHM13v2, GRCm39
from stories.RIP import annotation

SAVETO = ld.TABLES
SAVETO.mkdir(parents=True, exist_ok=True)

pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

for config in annotation.ld.Config.load():
    df = pd.read_pickle(config.summary)

    # Resolve the annotation based on RIPs from infected cells
    refexp = ld.c.refexp[config.ind]
    repcls = {"Z-RNA[CHM13v2]": CHM13v2.repcls, "Z-RNA[GRCm39]": GRCm39.repcls}[config.ind]
    df = ld.resolve.by_refexp(df, refexp, repcls)

    # Drop loop size for "unresolved" clusters
    df.loc[df['sequence[cls]'] == 'Unresolved', 'loop-size'] = -1

    # Use manual annotation where available and automatic annotation otherwise
    df['location'], df['Gene name'], df['Ensembl ID'] = \
        df['location [manual]'], df['Gene name [manual]'], df['Ensembl ID [manual]']

    df['location [auto]'] = df['location [auto]'].fillna('intergenic')

    mask = df['location'].isnull()
    df.loc[mask, 'location'] = df.loc[mask, 'location [auto]']
    df.loc[mask, 'Gene name'] = df.loc[mask, 'Gene name [auto]']
    df.loc[mask, 'Ensembl ID'] = df.loc[mask, 'Ensembl ID [auto]']

    # Mark clusters induced during the IAV/HSV-1 infection (in Z22 RIP or input)
    inducible = ld.c.inducible[config.ind]
    df = ld.resolve.virus_induction(
        df,
        inducible,
        log2fc=ld.c.thresholds.induced.log2fc,
        padj=ld.c.thresholds.induced.padj
    )

    # Mark clusters enriched in relevant RIPs
    enriched = ld.c.enriched[config.ind]
    df = ld.resolve.rip_enrichment(
        df, enriched, log2fc=ld.c.thresholds.enriched.log2fc, padj=ld.c.thresholds.enriched.padj
    )

    # Mark enriched & induced clusters
    for virus in "HSV-1", "IAV":
        for ab in "Z22", "FLAG":
            df[f"{virus} induced & {ab} enriched"] = df[f"{virus} induced"] & df[f"{ab} enriched [{virus}]"]

    # Simplify the location annotation
    df['location'] = df['location'].apply(lambda x: ld.c.locations.mapping.get(x, x))
    assert set(df['location']) <= set(ld.c.locations.palette), set(df['location']) - set(ld.c.locations.palette)

    # Mark dsRNA clusters as edited if they overlap with at least one edited sequence
    df['A-I edited'] = df['A2I'] > 0

    # Subsample the columns
    columns = list(chain(*[
        (name, (name, 'log2FoldChange'), (name, 'padj'), (name, 'pvalue')) for name in enriched
    ]))

    df = df[[
        'Ensembl ID', 'Gene name', 'location', 'loop-size', 'sequence[cls]', 'sequence[family]', 'A-I edited',
        'RIP-qPCR', 'HSV-1 induced', 'IAV induced', *columns,
        'HSV-1 induced & Z22 enriched', 'IAV induced & Z22 enriched',
        'HSV-1 induced & FLAG enriched', 'IAV induced & FLAG enriched',
    ]]

    # Save the tables
    df.to_pickle(SAVETO / f"{config.ind}.pkl")
    df.to_csv(SAVETO / f"{config.ind}.csv")
