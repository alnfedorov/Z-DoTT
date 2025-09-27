import pandas as pd


def resolve(df: pd.DataFrame, padj: float, log2fc: float, delta: float) -> pd.DataFrame:
    df['delta'] = df['Score [trt]'] - df['Score [ctrl]']

    df['category'] = 'Not significant'

    significant_up = (df['padj'] <= padj) & (df['log2fold_change'] >= log2fc)
    df.loc[significant_up, 'category'] = 'Significant up'

    strong_up = significant_up & (df['delta'] >= delta)
    df.loc[strong_up, 'category'] = 'Strong up'

    significant_down = (df['padj'] <= padj) & (df['log2fold_change'] <= -log2fc)
    df.loc[significant_down, 'category'] = 'Significant down'

    strong_down = significant_down & (df['delta'] <= -delta)
    df.loc[strong_down, 'category'] = 'Strong down'
    return df
