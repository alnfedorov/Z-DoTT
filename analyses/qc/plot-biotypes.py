import math

import matplotlib.pyplot as plt
import pandas as pd

import ld

ld.BIOTYPES.PLOTS.mkdir(parents=True, exist_ok=True)
plt.rcParams['svg.fonttype'] = 'none'


def pie_charts(counts: pd.DataFrame, palette: dict[str, str], order: list[str], postfix: str):
    for prj, df in counts.groupby('Project'):
        assert df['Title'].nunique() == len(df), f"Non-unique titles: {df['Title'].value_counts()}"

        melted = df.melt(id_vars=['Title'], value_vars=order, var_name='Category', value_name='Fragments')
        titles = df['Title'].unique()

        rows = int(math.ceil(len(titles) / 6))
        fig, axes = plt.subplots(rows, 6, figsize=(6 * 5, 5 * rows))
        axes = list(axes.ravel())

        melted = melted.sort_values(by='Title')
        for title, subdf in melted.groupby('Title'):
            ax = axes.pop()
            subdf['Color'] = subdf['Category'].apply(lambda x: palette[x])
            subdf['order'] = subdf['Category'].apply(lambda x: order.index(x))
            subdf = subdf.sort_values(by='order')

            ax.pie(
                subdf['Fragments'], pctdistance=0.75, labels=subdf['Category'], colors=subdf['Color'],
                autopct='%1.1f%%',
                textprops=dict(fontsize='small'), startangle=90,
                wedgeprops=dict(width=0.5, linewidth=0.9, edgecolor='black')
            )
            ax.text(0.5, 0.5, title, va='center', ha='center', fontsize='small', fontweight='bold',
                    transform=ax.transAxes)

            mlnfrag = subdf['Fragments'].sum() / 1_000_000
            ax.text(0.5, 0, f"N={mlnfrag:.1f}mln", va='top', ha='center', transform=ax.transAxes)
        ax.legend()

        fig.suptitle(prj, fontsize='xx-large', y=0.95, fontweight='bold')
        saveto = ld.BIOTYPES.PLOTS / f"{prj.replace('/', '-')}.{postfix}.svg"
        fig.savefig(saveto, bbox_inches="tight", pad_inches=0, transparent=True, dpi=300)
        # fig.show()
        plt.close(fig)


# Collapse basic categories
counts = pd.read_pickle(ld.BIOTYPES.COLLAPSED_FRAGMENTS)
mapping = {
    'mRNA': [
        '5\'UTR', 'CDS', '3\'UTR'
    ],
    'Other RNA': [
        'pseudogene', 'lncRNA', 'ribozyme', 'MT',
        'TEC', 'mRNA defective', 'rRNA', 'artifact', 'misc_RNA', 'processed_transcript'
    ],
    'smRNA': [
        'smRNA', '7SL', '7SK'
    ],
    'vRNA': [
        'HSV-1(+)', 'IAV(+)',
        'HSV-1(-)', 'IAV(-)',
    ]
}
for newcol, oldcols in mapping.items():
    oldcols = [col for col in oldcols if col in counts.columns]
    counts[newcol] = counts[oldcols].sum(axis=1)
    counts.drop(columns=oldcols, inplace=True)

# Reindex projects and experiments
experiments = {}
for series in ld.BIOTYPES.series.values():
    for prj in series:
        for exp in prj.experiments:
            experiments[(prj.ind, exp.ind)] = exp

counts['Experiment'] = [experiments[key] for key in counts['source']]

counts['Title'] = counts['Experiment'].apply(lambda x: x.attributes['title'] + '+' + x.ind)
counts['Treatment'] = counts['Experiment'].apply(lambda x: x.sample.attributes['treatment'])
counts['Replica'] = counts['Experiment'].apply(lambda x: x.sample.attributes['replica'])

selection = []
for exp in counts['Experiment']:
    if exp.library.selection == {'Total RNA', 'rRNA depletion'}:
        selection.append('input')
    else:
        s = exp.library.selection - {'Total RNA', 'rRNA depletion'}
        assert len(s) == 1
        selection.append(s.pop())
counts['Selection'] = selection
counts['Project'] = counts['source'].apply(lambda x: x[0])
counts = counts.drop(columns=['Experiment', 'source'])

order = [
    'vRNA', 'mRNA', 'intron', 'intergenic', 'Other RNA'
]
assert set(order) == set(counts.columns) - {'Assembly', 'Title', 'Treatment', 'Replica', 'Selection', 'Project'}

# Add pooled data
pooled = []
counts['group'] = counts['Title'].apply(
    lambda x: "_".join(x.split("_")[:-1])
    .replace("FLAG-input", "input").replace("Z22-input", "input")
)

for (_, title), subdf in counts.groupby(['Project', 'group']):
    record = {'Replica': 'pooled'}
    for col in ['Assembly', 'Project', 'Treatment', 'Selection']:
        assert subdf[col].nunique() == 1, f"Non-unique {col}: {subdf[col].unique()}"
        record[col] = subdf[col].iloc[0]

    for col in order:
        record[col] = subdf[col].mean()
    record['Title'] = f"{title}_pooled"
    pooled.append(record)
counts = pd.concat([counts, pd.DataFrame(pooled)]).drop(columns='group')

palette = dict(zip(order, plt.colormaps["tab10"].colors))
pie_charts(counts, palette, order, "biotypes")

# Plot only cellular and vRNAs
order.remove('vRNA')
counts['Cellular RNA'] = counts[order].sum(axis=1)
counts.drop(columns=order, inplace=True)
order = ['vRNA', 'Cellular RNA']
palette = {
    "vRNA": "#99D4C0",
    "Cellular RNA": "#FBB273"
}
pie_charts(counts, palette, order, "cellular-vRNA")
