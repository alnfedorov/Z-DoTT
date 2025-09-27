import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import ld
from stories import normalization

plt.rcParams['svg.fonttype'] = 'none'
ld.PCA.mkdir(parents=True, exist_ok=True)

counts = pd.read_pickle(normalization.FRAGMENTS)

# Parse titles for all experiments
titles = {}
for projects in normalization.SERIES.values():
    for prj in projects:
        for exp in prj.experiments:
            titles[(prj.ind, exp.ind)] = exp.attributes['title']

for cnts in counts.values():
    cnts['Title'] = cnts['source'].map(titles)

# Re-group by project
groups = {}
for cnts in counts.values():
    prjs = cnts['source'].apply(lambda x: x[0]).unique()
    for prj in prjs:
        subdf = cnts[cnts['source'].apply(lambda x: x[0]) == prj]
        assert prj not in groups
        groups[prj] = subdf.copy()

# Extract extra fields from the title for selected projects
for prjind, mappings in [
    ('MEF IAV batch 1 [B256178, JCC400]', {"condition": 1, "IP": 2}),
    ('MEF HSV-1 batch 1 [B261790, JCC411]', {"condition": 1, "IP": 2}),
    ('MEF HSV-1 batch 2 [B319096, JCC478]', {"condition": 1, "IP": 2}),
    ('HT-29 HSV-1/IAV [B831009, JCC411]', {"condition": 1, "IP": 2}),
]:
    if prjind not in groups:
        print(f"Skipping {prjind}")
        continue

    subdf = groups[prjind]
    for field, index in mappings.items():
        subdf[field] = subdf['Title'].str.split('_').apply(
            lambda x: x[index].replace("FLAG-input", "input").replace("Z22-input", "input")
        )


def process(title: str, subdf: pd.DataFrame, hue: str | None = None, style: str | None = None):
    # Ignore zeros
    columns = (subdf.select_dtypes('number') > 1e-12).any()
    columns = columns[columns].index.tolist()

    # Normalize for the sequencing depth
    total = subdf[columns].sum(axis=1)
    subdf.loc[:, columns] = subdf[columns].div(total, axis=0)

    # Run PCA
    pipeline = Pipeline([('scaling', StandardScaler()), ('pca', PCA(n_components=2))])
    transformed = pipeline.fit_transform(subdf[columns]).T
    subdf['PCA-1'] = transformed[0]
    subdf['PCA-2'] = transformed[1]

    # Make the plot
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.scatterplot(
        data=subdf, x='PCA-1', y='PCA-2', hue=hue, style=style, ax=ax,
        palette={"Z22-RIP": "#1f77b4", "FLAG-RIP": "#ff7f0e", "input": "#2ca02c"}
    )

    for smpl, (x, y) in zip(subdf['Title'], zip(transformed[0], transformed[1])):
        ax.text(x, y, smpl)

    ax.set_title(title)
    # ax.grid()
    ax.spines[['top', 'right']].set_visible(False)

    explained = pipeline['pca'].explained_variance_ratio_
    ax.set_xlabel(f"PCA-1 ({explained[0]:.1%})")
    ax.set_ylabel(f"PCA-2 ({explained[1]:.1%})")

    title = title.replace(" ", "-").replace("/", "-")
    saveto = ld.PCA / f"{title}.svg"
    fig.savefig(saveto, bbox_inches="tight", pad_inches=0, transparent=True)
    fig.show()
    plt.close(fig)


# Pseudo project for all MEFs datasets
pooled = []
for prefix, prjind in [
    ("B1", "MEF IAV batch 1 [B256178, JCC400]"),
    ("B1", "MEF HSV-1 batch 1 [B261790, JCC411]"),
    ("B2", "MEF HSV-1 batch 2 [B319096, JCC478]"),
]:
    if prjind not in groups:
        print(f"Skipping {prjind}")
        continue
    subdf = groups[prjind].copy()
    subdf['Title'] = subdf['Title'].apply(lambda x: f"{prefix}: {x}")
    pooled.append(subdf)
groups['MEF'] = pd.concat(pooled)

# Drop all IgG RIPs
for prjind, df in groups.items():
    if 'IP' in df:
        groups[prjind] = df[~df['IP'].isin({'IgG', 'IgG-RIP'})].copy()

# PCA for all data
for prjind, title, hue, style in [
    ('MEF', "MEF: All", "IP", "condition"),
    ('MEF IAV batch 1 [B256178, JCC400]', "MEF: IAV [N1]", "IP", "condition"),
    ('MEF HSV-1 batch 1 [B261790, JCC411]', "MEF: HSV-1 [N1]", "IP", "condition"),
    ('MEF HSV-1 batch 2 [B319096, JCC478]', "MEF: HSV-1 [N2]", "IP", "condition"),
    ('HT-29 HSV-1/IAV [B831009, JCC411]', "HT-29: IAV+HSV-1 [N1]", "IP", "condition"),
]:
    if prjind not in groups:
        print(f"Skipping {prjind}")
        continue
    process(title, groups[prjind], hue, style)
