import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.stats import pearsonr

import ld
from stories.RIP.plots.ld import TABLES as ZRNA_TABLES

PALETTE = {
    'None / Weak': '#99dbff',
    'HSV-1 specific': '#ff3333',
    'IAV specific': '#5cad5c',
    'HSV-1 and IAV': '#c58544'
}

ALPHA = {
    'None / Weak': 0.25,
    'HSV-1 specific': 1.0,
    'IAV specific': 1.0,
    'HSV-1 and IAV': 1.0
}

SIZE = {
    'None / Weak': 5,
    'HSV-1 specific': 10,
    'IAV specific': 10,
    'HSV-1 and IAV': 10
}

workloads = pd.read_pickle(ld.WORKLOAD)
for name, iav, hsv in [
    ("Z-RNA[GRCm39]", "MEF: IAV vs mock", "MEF: HSV-1 vs mock"),
    ("Z-RNA[CHM13v2]", "HT-29: IAV vs mock", "HT-29: HSV-1 vs mock")
]:
    # Select and merge estimated downstream transcription rates for IAV and HSV-1 infections
    iav = workloads[iav].aberrantome["Downstream"][['Ensembl ID', 'Name', 'Score [trt]', 'label']]
    hsv = workloads[hsv].aberrantome["Downstream"][['Ensembl ID', 'Name', 'Score [trt]', 'label']]

    df = pd.merge(iav, hsv, on=['Ensembl ID', 'Name'], suffixes=("_IAV", "_HSV"), how='inner')
    df = df.rename(columns={
        'Score [trt]_IAV': 'IAV Downstream Transcript Rate',
        'Score [trt]_HSV': 'HSV-1 Downstream Transcription Rate'
    })

    assert (df['label_HSV'].isna() | (df['label_HSV'] == df['label_IAV'])).all()
    df['label'] = df['label_IAV']
    df = df.drop(columns=['label_IAV', 'label_HSV'])

    # Link them to detected Z-RNAs in HSV-1/IAV infected cells
    zrna = pd.read_pickle(ZRNA_TABLES / f"{name}.pkl")

    for category, z22, flag in [
        ("Z-RNA[HSV-1]", "HSV-1 induced & Z22 enriched", "HSV-1 induced & FLAG enriched"),
        ("Z-RNA[IAV]", "IAV induced & Z22 enriched", "IAV induced & FLAG enriched")
    ]:
        passed = set()
        for column in z22, flag:
            mask = zrna[column] & zrna['Ensembl ID'].notna() & (zrna['location'] == 'Elongated 3` end')
            print(f'\t{column} with N={mask.sum()} Z-RNAs in Elongated 3` end')
            passed |= set(zrna.loc[mask, 'Ensembl ID'])
        passed = {x.split('.')[0] for x in passed}

        df[category] = df['Ensembl ID'].map(lambda x: x.split('.')[0] in passed)
        print(f"\t{category} -> {df[category].sum()}\n")

    # Categorize all hits based on the Z-RNAs downstream
    df['category'] = 'None / Weak'
    df.loc[df['Z-RNA[HSV-1]'], 'category'] = 'HSV-1 specific'
    df.loc[df['Z-RNA[IAV]'], 'category'] = 'IAV specific'
    df.loc[df['Z-RNA[HSV-1]'] & df['Z-RNA[IAV]'], 'category'] = 'HSV-1 and IAV'

    # Make the plot
    X = 'HSV-1 Downstream Transcription Rate'
    Y = 'IAV Downstream Transcript Rate'
    limit = 0.5

    jointplot = sns.JointGrid(
        data=df, x=X, y=Y, xlim=(0, 1), ylim=(0, 1), ratio=13, hue='category', palette=PALETTE
    )

    # Marginal plots
    jointplot.plot_marginals(
        sns.histplot, common_norm=False, common_bins=True, binwidth=0.01, binrange=(0, 1),
        stat="probability", fill=False, kde=True, alpha=0.0,
        kde_kws={"gridsize": 500, "clip": (0, 1)}
    )
    jointplot.ax_marg_x.set(xlim=(0, limit), ylim=(0, 0.1))
    jointplot.ax_marg_y.set(ylim=(0, limit), xlim=(0, 0.1))

    # Joint plot
    for category, subdf in df.groupby('category'):
        jointplot.ax_joint.scatter(
            subdf[X], subdf[Y],
            alpha=ALPHA[category], s=SIZE[category], linewidth=0, color=PALETTE[category]
        )

    # Regression line
    sns.regplot(
        data=df, x=X, y=Y, ax=jointplot.ax_joint, truncate=False,
        scatter=False, color='black', ci=95, line_kws={"lw": 1}, robust=True
    )

    offset = ld.config.threshold.delta
    jointplot.ax_joint.plot([offset, 1.0], [0.0, 1 - offset], color='black', ls='--', lw=1.5)
    jointplot.ax_joint.plot([0.0, 1 - offset], [offset, 1.0], color='black', ls='--', lw=1.5)

    points = ((df[X] - df[Y]).abs() <= offset).mean()
    jointplot.ax_joint.text(
        0.99, 0.99, f'|diff| <= {offset} = {points:.1%}',
        va='top', ha='right', transform=jointplot.ax_joint.transAxes
    )

    # Ticks and limits
    jointplot.ax_joint.set(
        xlim=(0, limit), ylim=(0, limit),
        xticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5], xticklabels=[0, 10, 20, 30, 40, 50],
        yticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5], yticklabels=[0, 10, 20, 30, 40, 50],
    )

    # Legend
    cnts = df['category'].value_counts()
    legend = []
    for key, color in sorted(PALETTE.items()):
        legend.append(Line2D(
            [0], [0], marker='o', color='white', label=f"[N={cnts.get(key, 0)}] {key}",
            markerfacecolor=color, markersize=8
        ))
    jointplot.ax_joint.legend(handles=legend, loc='upper left', labelspacing=0.1, handletextpad=0.02, alignment='left')

    # Labels
    # points = ((df[X] - df[Y]).abs() <= offset).mean()
    # jointplot.ax_joint.text(1.0, 0.04, f'|diff| <= {offset} = {points:.1%}', va='bottom', ha='right')

    pearson = pearsonr(df[X], df[Y])
    pearson = pearson.statistic
    jointplot.ax_joint.text(
        1.0, 0.02, f'r={pearson:.3f}', va='bottom', ha='right',
        transform=jointplot.ax_joint.transAxes
    )

    # Title
    jointplot.fig.suptitle(name)

    # Annotate relevant targets
    labeled = df['label'].notna()
    print(f"Annotating {labeled.sum()} hits")
    for _, hit in df[labeled].iterrows():
        jointplot.ax_joint.text(
            hit[X], hit[Y], hit['label'], fontsize=12, ha='left', va='bottom', fontdict={'weight': 'bold'}
        )
        print(f"\t{hit['label']} -> {hit[X]:.2f} vs {hit[Y]:.2f} [{hit['category']}]")
    print()

    # jointplot.fig.show()
    for format in ['png', 'svg']:
        saveto = ld.RESULTS / f"Z-RNA concordance {name}.{format}"
        jointplot.fig.savefig(saveto, bbox_inches="tight", pad_inches=0, transparent=True, dpi=400)
    plt.close(jointplot.fig)
