from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from biobit.toolkit import seqproj
from matplotlib import transforms

import ld
from stories.nextflow import series

plt.rcParams['svg.fonttype'] = 'none'

pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

ld.c.plots.saveto.mkdir(parents=True, exist_ok=True)


def bar_plot(
        df: pd.DataFrame, *, x: str, order: list[str], hue: str, hue_order: list[str], palette: dict[str, str],
        ylim: int = 100, title: str, saveto: Path
):
    order = [i for i in order if i in df[x].unique()]
    hue_order = [i for i in hue_order if i in df[hue].unique()]

    scale = 7
    fig, ax = plt.subplots(figsize=((20 + 20 * len(order)) / scale, 150 / scale))

    sns.stripplot(
        data=df, x=x, y='sites', hue=hue, ax=ax, order=order, palette=palette,
        hue_order=hue_order, dodge=True, size=27.5, alpha=1.0, linewidth=1.25, jitter=False
    )
    sns.barplot(
        data=df, x=x, y='sites', hue=hue, ax=ax, order=order, lw=2.0, palette=palette,
        hue_order=hue_order, errorbar=("sd", 1), capsize=0.3, edgecolor="black",
        err_kws={'linewidth': 3.5, 'color': 'black'}
    )

    # Annotate the total number of sites
    total = df.groupby(['x', 'ind'], as_index=False)['total'].sum() \
        .groupby('x')['total'].median().astype(int).to_dict()
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    for x, v in enumerate(order):
        ax.text(x, 1.0, f"{total[v]:,}", transform=trans, size=5 * scale, va='bottom', ha='center')
    ax.text(-0.045, 1.0, 'N:', transform=ax.transAxes, size=5 * scale, va='bottom', ha='right')

    ax.set_title(title, y=1.05, fontsize=6 * scale)
    ax.set_xlabel('')
    ax.set_ylabel('Localization (%)', fontsize=6 * scale)
    ax.tick_params(axis='both', which='major', labelsize=5 * scale, width=0.35 * scale, length=1.85 * scale)
    ax.set_yticks(range(0, int(ylim) + 1, 5))
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=5 * scale)

    ax.spines[['right', 'top']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(0.35 * scale)
    ax.set_ylim(0, ylim)

    fig.show()
    fig.savefig(saveto, bbox_inches="tight", transparent=True, pad_inches=0)
    plt.close(fig)


# RIP data #############################################################################################################
for title, projects, ylim in [
    ("MEF: IAV", [series.internal.B256178], 27.5),
    ("HT-29: IAV", [series.internal.B831009], 27.5),
    ("MEF: HSV-1", [series.internal.B261790, series.internal.B319096], 32.5),
    ("HT-29: HSV-1", [series.internal.B831009], 32.5),
]:
    def meta(exp: seqproj.Experiment) -> dict[str, str]:
        selection = exp.attributes['title'].split('_')[-2].replace("-RIP", "")
        if "input" in selection:
            selection = "input"
        return {
            "treatment": exp.sample.attributes["treatment"], "replica": exp.sample.attributes["replica"],
            "ind": exp.ind, "selection": selection
        }


    assembly = "GRCm39" if "Mus musculus" in projects[0].samples[0].organism else "CHM13v2"
    experiments = {prj.ind: prj.experiments for prj in projects}
    df = ld.utils.load.annotated(
        ld.c.reat.annotated / assembly, experiments, meta, ld.c.plots.mapping,
        pivot=["treatment", "replica", "ind", "selection"]
    )
    df['x'] = [f"{treatment} {selection}" for treatment, selection in zip(df['treatment'], df['selection'])]
    order = [
        f"{treatment} {selection}"
        for treatment in ["mock", "HSV-1", "IAV"]
        for selection in ["input", "FLAG", "Z22"]
        if treatment in df['treatment'].unique() and selection in df['selection'].unique()
    ]

    prjind = "+".join([prj.ind.replace(" ", "_").replace("/", "_") for prj in projects])
    bar_plot(
        df, x='x', order=order, palette=ld.c.plots.palette,
        hue='location', hue_order=["Proximal intergenic", "Distal intergenic"], ylim=ylim,
        saveto=ld.c.plots.saveto / f"{title}.svg", title=f"{title}\n{prjind}"
    )
