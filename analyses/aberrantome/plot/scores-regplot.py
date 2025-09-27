from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from joblib import Parallel, delayed
from matplotlib.lines import Line2D
from scipy.stats import pearsonr

import ld

ALPHA = {
    "Not significant": 0.4,
    "Significant up": 0.6,
    "Strong up": 0.6,
    "Significant down": 0.6,
    "Strong down": 0.6
}


def job(df: pd.DataFrame, title: str, score: str, trtlabel: str, ctrlabel: str, saveto: Path):
    plt.rcParams['svg.fonttype'] = 'none'
    X, Y = f"Score [trt]", f"Score [ctrl]"
    limit = ld.c.threshold.trlimit

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.spines[['top', 'right']].set_visible(False)

    # Scatter plot
    for category, subdf in df.groupby('category'):
        ax.scatter(
            subdf[X], subdf[Y], s=11, edgecolor='black', alpha=ALPHA[category],
            linewidth=0.25, color=ld.c.palette[category]
        )

    # Ticks and limits
    ticks = [0.1 * i for i in range(10) if 0.1 * i <= limit]
    labels = [f"{int(100 * i)}" for i in ticks]
    ax.set(
        xlim=(0, limit), xticks=ticks, xticklabels=labels,
        ylim=(0, limit), yticks=ticks, yticklabels=labels,
    )

    # Regression line
    sns.regplot(
        data=df, x=X, y=Y, ax=ax, truncate=False,
        scatter=False, color='black', ci=95, line_kws={"lw": 1}, robust=False
    )

    # X/Y Labels
    ax.set(
        xlabel=f"{trtlabel} {score} transcription rate (%)",
        ylabel=f"{ctrlabel} {score} transcription rate (%)"
    )

    # 'Strong' changes thresholds
    offset = ld.config.threshold.delta
    ax.plot([offset, limit], [0.0, limit - offset], color='black', ls='--', lw=1.5)
    ax.plot([0.0, limit - offset], [offset, limit], color='black', ls='--', lw=1.5)

    # Legend
    cnts = df['category'].value_counts()
    for key, loc in [
        ("Significant up", [0.75, 0.3]), ("Significant down", [0.3, 0.75]),
        ("Strong up", [0.75, 0.25]), ("Strong down", [0.3, 0.7])
    ]:
        value, color = cnts.get(key, 0), ld.c.palette[key]
        ax.text(
            loc[0], loc[1], f"N={value}", va='bottom', ha='left', fontsize=12, color=color, weight='bold',
            transform=ax.transAxes
        )

    legend = []
    for key, color in ld.c.palette.items():
        legend.append(Line2D(
            [0], [0], marker='o', color='white', label=f"[N={cnts.get(key, 0)}] {key}",
            markerfacecolor=color, markersize=8
        ))
    ax.legend(handles=legend, loc='upper left', labelspacing=0.1, handletextpad=0.02, alignment='left')

    # Stats
    points = ((df[X] - df[Y]).abs() <= offset).mean()
    ax.text(1.0, 1.0, f'|diff| <= {offset} = {points:.1%}', va='top', ha='right', transform=ax.transAxes)

    pearson = pearsonr(df[X], df[Y])
    pearson = pearson.statistic
    ax.text(
        1.0, 0.02, f'r={pearson:.3f}', va='bottom', ha='right',
        transform=ax.transAxes
    )

    # Title
    fig.suptitle(title)

    # # Annotate top hits
    # sig['order'] = (sig[X] - sig[Y]).abs()
    # sig = sig.sort_values('order', ascending=False)
    # for top in sig[sig[X] - sig[Y] > 0], sig[sig[X] - sig[Y] < 0]:
    #     top = top.head(45)
    #     for x, y, gene in zip(top[X], top[Y], top['Name']):
    #         ax.text(x, y, gene, fontsize=7, ha='center', va='bottom')

    # Annotate relevant targets
    labeled = df['label'].notna()
    for _, hit in df[labeled].iterrows():
        match hit['category']:
            case 'Strong up' | 'Strong down':
                ax.text(
                    hit[X], hit[Y], hit['label'], fontsize=12, ha='left', va='bottom', fontdict={'weight': 'bold'}
                )
            case 'Significant up' | 'Significant down':
                ax.text(
                    hit[X], hit[Y], hit['label'], fontsize=12, ha='left', va='bottom', fontdict={'weight': 'normal'}
                )
            case _:
                ax.text(
                    hit[X], hit[Y], hit['label'], fontsize=7, ha='left', va='bottom', fontdict={'weight': 'light'}
                )

    # Save the result
    saveto.parent.mkdir(parents=True, exist_ok=True)
    for st in saveto.with_suffix(".png"), saveto.with_suffix(".svg"):
        fig.savefig(st, bbox_inches="tight", pad_inches=0, transparent=True, dpi=400)

    # Close the figure
    plt.close(fig)


workloads = []
for ind, wd in pd.read_pickle(ld.WORKLOAD).items():  # type: str, ld.config.Workload
    for name, df in wd.aberrantome.items():
        title = f"{ind} ({name})"
        saveto = ld.RESULTS / ind / f"{name} regplot"
        workloads.append((df, title, name, wd.comparison.treatment_label, wd.comparison.control_label, saveto))

Parallel(n_jobs=-1)(delayed(job)(*args) for args in workloads)
