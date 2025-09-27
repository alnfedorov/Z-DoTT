import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import pandas as pd
import seaborn as sns
from biobit.toolkit import seqproj

import ld
from stories.nextflow import series

plt.rcParams['svg.fonttype'] = 'none'

pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

ld.c.plots.saveto.mkdir(parents=True, exist_ok=True)


def series_plot(df: pd.DataFrame, *, x: str, order: list[str], saveto: Path, ylim: int = 100):
    order = [item for item in order if item in set(df[x])]
    df = df[df[x].isin(order)].copy()
    df['_x'] = df[x].apply(lambda x: order.index(x))

    scale = 7
    fig, ax = plt.subplots(figsize=(180 / scale, 106 / scale))

    ax = sns.scatterplot(data=df, x='_x', y='sites', color='black', legend=False, s=90 * scale, zorder=1000, ax=ax)
    ax = sns.lineplot(
        data=df, x='_x', hue='location', y='sites', lw=1 * scale,
        hue_order=["Exonic", "Intronic", "Proximal intergenic", "Distal intergenic"],
        palette=ld.c.plots.palette, ax=ax
    )

    ax.set_ylabel("Localization (%)", fontsize=6 * scale)
    ax.set_xlabel(f"Time point", fontsize=6 * scale)

    ax.set_xticks(range(len(order)), order, fontsize=5 * scale)  # , rotation=15)
    ax.tick_params(axis='both', which='major', labelsize=5 * scale, width=0.35 * scale, length=1.85 * scale)
    ax.spines[['left', 'bottom']].set_linewidth(0.35 * scale)
    ax.spines[['right', 'top']].set_visible(False)

    # Annotate the total number of sites
    total = df.groupby([x])['total'].sum().to_dict()
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    for x, v in enumerate(order):
        ax.text(x, 1.0, f"{total[v]:,}", transform=trans, size=4 * scale, va='bottom', ha='center')
    ax.text(-0.045, 1.0, 'N:', transform=ax.transAxes, size=4 * scale, va='bottom', ha='right')

    # ax.grid(linewidth=0.35 * scale)

    # ax.set_title("A-to-I editing sites distribution", fontsize=6 * scale)

    ax.set(ylim=(0, ylim), xlim=(0, len(order) - 1))
    legend = ax.legend(labelspacing=0, prop={'size': 6 * scale}, loc='upper right')
    for l in legend.legend_handles:
        l.set_linewidth(1.0 * scale)

    fig.savefig(saveto, bbox_inches="tight", transparent=True, pad_inches=0)
    fig.show()
    plt.close(fig)


# Load pooled experiments
pooled = {}
with open(ld.c.reat.seqproj, 'rb') as stream:
    for prjind, experiments in pickle.load(stream).items():
        pooled[prjind] = [exp for exp in experiments if exp.sample.attributes['replica'] == 'pooled']


# HSV-1 time-course ####################################################################################################
def meta(exp: seqproj.Experiment) -> dict[str, str]:
    return {
        "time-point": exp.sample.attributes['time-point'], "protocol": exp.sample.attributes["protocol"],
        "replica": exp.sample.attributes["replica"]
    }


project = series.SRA.PRJNA256013
df = ld.utils.load.annotated(
    ld.c.reat.annotated / "CHM13v2", {project.ind: pooled[project.ind]}, meta, ld.c.plots.mapping,
    pivot=["replica", "time-point", "protocol"]
)

prjind = project.ind.replace(" ", "_").replace("/", "_")
# Time labels
order = [
    'mock', '0-1hpi', '1-2hpi', '2-3hpi', '2hpi', '3-4hpi', '4-5hpi', '4hpi',
    '5-6hpi', '6-7hpi', '6hpi', '7-8hpi', '8hpi',
]
assert set(order) == set(df['time-point'])

subdf = df[(df['protocol'] == '4sU-RNA') & (df['replica'] == "pooled")].copy()
series_plot(subdf, x='time-point', order=order, saveto=ld.c.plots.saveto / "HSV-1 time course(4sU).svg", ylim=65)

subdf = df[(df['protocol'] == 'Total-RNA') & (df['replica'] == 'pooled')].copy()
series_plot(subdf, x='time-point', order=order, saveto=ld.c.plots.saveto / "HSV-1 time course(Total).svg", ylim=60)


# IAV time-course ######################################################################################################
def meta(exp: seqproj.Experiment) -> dict[str, str]:
    return {
        "cells": exp.sample.attributes['cells'], "IAV": exp.sample.attributes['IAV'],
        "time-point": exp.sample.attributes['time-point'], "replica": exp.sample.attributes["replica"]
    }


project = series.SRA.PRJNA382632
df = ld.utils.load.annotated(
    ld.c.reat.annotated / "CHM13v2", {project.ind: pooled[project.ind]}, meta, ld.c.plots.mapping,
    pivot=["cells", "IAV", "time-point", "replica"]
)

# Average different mock time points
mock = df[df['IAV'] == 'mock'].copy()
mock['time-point'] = 'mock'

index = ['cells', 'IAV', 'time-point', 'replica', 'location']
mock = mock.groupby(index)['total'].median().reset_index()
total = mock['total'].sum()
mock['sites'] = mock['total'] / total * 100


df = pd.concat([mock, df[df['IAV'] != 'mock']])

for IAV in ['H1N1', 'H3N2', 'H5N1']:
    subdf = df[df['IAV'].isin({IAV, 'mock'})].copy()
    order = [x for x in ('mock', '3h', '6h', '12h') if x in subdf['time-point'].unique()]

    series_plot(
        subdf, x='time-point', order=order, saveto=ld.c.plots.saveto / f"IAV time course(MDM + {IAV}).svg", ylim=47
    )
