import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import ld
from utils.plot import VolcanoPreset

plt.rcParams['svg.fonttype'] = 'none'

for table in ld.TABLES.glob("*.pkl"):
    ind = table.stem
    df = pd.read_pickle(table)

    saveto = ld.RESULTS / ind
    saveto.mkdir(exist_ok=True, parents=True)

    for virus in "HSV-1", "IAV":
        # Select only induced hits
        induced = df[df[f'{virus} induced']].copy()
        for ab in "Z22", "FLAG":
            title = f"{ab} enriched [{virus}]"

            # induced.loc[induced[(title, 'pvalue')] == 0, (title, 'pvalue')] = \
            #     induced.loc[induced[(title, 'pvalue')] != 0, (title, 'pvalue')].min()

            induced['-log10pv'] = -np.log10(induced[(title, 'pvalue')])
            induced['log2FoldChange'] = induced[(title, 'log2FoldChange')]
            induced['padj'] = induced[(title, 'padj')]

            induced['category'] = induced['location']
            induced.loc[~induced[title], 'category'] = 'Not significant'

            fig, ax = plt.subplots(figsize=(6, 6))

            xmax, ymax = ld.c.volcano.xymax[ind][title]
            VolcanoPreset(induced, '-log10pv', 'padj', 'log2FoldChange') \
                .threshold(padj=ld.c.thresholds.enriched.padj, log2fc=ld.c.thresholds.enriched.log2fc) \
                .color('category', ld.c.locations.palette) \
                .edge_width('A-I edited', {True: 0.5, False: 0}) \
                .size('A-I edited', {True: 13, False: 13}) \
                .xlimit(0, xmax) \
                .ylimit(0, ymax) \
                .set(alpha=0.55, annotate='RIP-qPCR') \
                .plot(ax)

            ax.set_title(f"{ind}: {title}", fontsize=12, loc='left', fontdict={'fontweight': 'bold'})
            ax.spines[['right', 'top']].set_visible(False)

            xlimits = ax.get_xlim()
            ax.set_xticks(range(int(xlimits[0]), int(xlimits[1]) + 1))

            # Make the legend
            legend, counts = [], induced['category'].value_counts().to_dict()
            for description in ld.c.locations.order:
                if counts.get(description, 0) > 0:
                    color = ld.c.locations.palette[description]
                    legend.append(plt.Line2D(
                        [0], [0], marker='o', color='white', markerfacecolor=color, markersize=10,
                        label=f"[N={counts[description]}] {description}",
                    ))

            edited = induced[title] & induced['A-I edited']
            edtotal, edfreq = edited.sum(), edited.sum() / induced[title].sum()
            legend.append(plt.Line2D(
                [0], [0], marker='o', color='black', label=f"A->I edited: {edtotal} ({edfreq:.1%})",
                markerfacecolor='white', markersize=10, markeredgewidth=0.5
            ))
            ax.legend(handles=legend, loc='upper left', title="Localization", labelspacing=0.1, handletextpad=0.02)

            fig.savefig(saveto / f"{title} volcano-plot.svg", bbox_inches="tight", pad_inches=0, transparent=True)
            # fig.show()
            plt.close(fig)
