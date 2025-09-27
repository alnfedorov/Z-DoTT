from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from matplotlib.lines import Line2D

import ld
from utils.plot import VolcanoPreset


def job(
        df: pd.DataFrame, title: str, xlabel: str,
        xlim: tuple[float | None, float | None], ylim: tuple[float | None, float | None],
        saveto: Path
):
    plt.rcParams['svg.fonttype'] = 'none'

    df = df.copy()
    df['category'] = df['category'].replace({
        "Strong up": "Significant up", "Strong down": "Significant down"
    })
    assert set(df['category']) <= {"Not significant", "Significant up", "Significant down"}

    # Replace 0 p-values with the smallest non-zero value and calculate -log10(pvalue)
    zeropval = df['pvalue'] == 0
    df.loc[zeropval, 'pvalue'] = df.loc[~zeropval, 'pvalue'].min()
    df['-log10(pvalue)'] = -np.log10(df['pvalue'])

    fig, ax = plt.subplots(figsize=(8, 7))

    VolcanoPreset(df, '-log10(pvalue)', 'padj', 'log2fc') \
        .threshold(padj=ld.c.threshold.padj, log2fc=ld.c.threshold.log2fc) \
        .color('category', ld.c.palette) \
        .ylimit(*ylim) \
        .xlimit(*xlim) \
        .set(alpha=0.55, annotate='label') \
        .plot(ax)

    ax.set_title(title, fontsize=16, loc='left', fontdict={'fontweight': 'bold'})
    ax.spines[['right', 'top']].set_visible(False)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("$log_{10}$(p-value)")
    xlimits = ax.get_xlim()
    ax.set_xticks(range(int(xlimits[0]), int(xlimits[1]) + 1))

    # Legend
    legend, counts = [], df['category'].value_counts().to_dict()
    for description, color in ld.c.palette.items():
        if counts.get(description, 0) > 0:
            legend.append(Line2D(
                [0], [0], marker='o', color='white', label=f"[N={counts[description]}] {description}",
                markerfacecolor=color, markersize=10,
            ))

    ax.legend(
        handles=legend, bbox_to_anchor=(-0.02, 1.01), loc='upper left',
        labelspacing=0.15, handletextpad=0.02, fontsize=12, frameon=False
    )

    # Save the figure
    saveto.parent.mkdir(parents=True, exist_ok=True)
    for st in saveto.with_suffix(".png"), saveto.with_suffix(".svg"):
        st.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(st, bbox_inches="tight", pad_inches=0, transparent=True, dpi=400)

    saveto.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(saveto, bbox_inches="tight", pad_inches=0, transparent=True, dpi=400)
    plt.close(fig)


workloads = []
for ind, wd in pd.read_pickle(ld.WORKLOAD).items():  # type: str, ld.config.Workload
    for name, df in wd.aberrantome.items():
        xlabel = f"{name} transcription rate"
        title = f"{ind} ({name})"
        saveto = ld.RESULTS / ind / f"{name} volcano"
        workloads.append((df, title, "$log_{2}$" + f"({xlabel})", (-8, 8), (0, 300), saveto))

Parallel(n_jobs=-1)(delayed(job)(*args) for args in workloads)
