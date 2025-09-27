import pickle

import pandas as pd
from joblib import delayed, Parallel

from stories.RIP import annotation


def make_summary(config: annotation.Config):
    partition, coordinates = [], []
    for element in config.elements:
        partition.append(element.ind)
        coordinates.append(f"{element.contig}:{element.envelope().start}-{element.envelope().end}")

    df = pd.DataFrame({"partition": partition, "coordinates": coordinates})

    with open(config.sequences, 'rb') as stream:
        mapping = pickle.load(stream)
    for col in [x.ind for x in config.comparisons]:
        df[f"sequence [{col}]"] = [mapping[col][x] for x in df["partition"]]

    with open(config.localization, 'rb') as stream:
        mapping = pickle.load(stream)
    for col in [x.ind for x in config.comparisons]:
        df[f"location-auto [{col}]"] = [mapping[col][x, "automatic"] for x in df["partition"]]
        df[f"location-manual [{col}]"] = [mapping[col][x, "manual"] for x in df["partition"]]

    with open(config.loops, 'rb') as stream:
        signal = pickle.load(stream)
    for col in [x.ind for x in config.comparisons]:
        df[f"loop-size [{col}]"] = [signal[col][x] for x in df["partition"]]

    for name, table in {
        "A2I": config.a2i,
        "RIP-qPCR": config.probes
    }.items():
        with open(table, 'rb') as stream:
            mapping = pickle.load(stream)
        df[name] = [mapping[x] for x in df["partition"]]

    with open(config.alltests, 'rb') as stream:
        test_results = pickle.load(stream)

    for name, table in test_results.items():
        for col in "log2FoldChange", "pvalue", "padj":
            mapping = table[col].to_dict()
            df[name, col] = [mapping.get(x, None) for x in df["partition"]]
        for col in "pvalue", "padj":
            df[name, col] = df[name, col].fillna(1.0)

    df = df.set_index(['partition', 'coordinates']).copy()
    df.to_pickle(config.summary)


Parallel(n_jobs=-1, backend='sequential')(delayed(make_summary)(config) for config in annotation.Config.load())
