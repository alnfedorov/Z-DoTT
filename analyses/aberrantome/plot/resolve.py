import pandas as pd

import ld
from stories.aberrantome.calculate import ld as calc
from stories.aberrantome.plot.ld.config import Workload

workload = {}
for cmp in calc.comparisons.load():
    assert cmp.ind not in workload, (cmp.ind, workload)

    aberrantome = {}
    for type, path in [
        ("Downstream", cmp.read_through), ("Upstream", cmp.read_in),
        ("Intronic", cmp.intronic), ("Divergent", cmp.divergent)
    ]:
        df = pd.read_csv(path)

        # Sanity check
        for col in ["Score [trt]", "Score [ctrl]"]:
            assert df[col].max() <= 1.0 and df[col].min() >= 0.0, df[col].describe()

        # Resolve stat. test categories
        df = ld.resolve(df, padj=ld.c.threshold.padj, log2fc=ld.c.threshold.log2fc, delta=ld.c.threshold.delta)

        # Annotate relevant genes
        df['label'] = pd.NA
        mask = df['Name'].isin(ld.config.annotate.by_assembly[cmp.assembly])
        df.loc[mask, 'label'] = df.loc[mask, 'Name']

        # Select the relevant columns
        aberrantome[type] = df[[
            'Ensembl ID', 'Name', 'Score [trt]', 'Score [ctrl]', 'pvalue', 'log2fold_change', 'padj',
            'category', 'label'
        ]].rename(columns={'log2fold_change': 'log2fc'})
    workload[cmp.ind] = Workload(cmp, aberrantome)

ld.WORKLOAD.parent.mkdir(parents=True, exist_ok=True)
pd.to_pickle(workload, ld.WORKLOAD)
