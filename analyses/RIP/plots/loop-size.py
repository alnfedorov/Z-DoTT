import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import ld

plt.rcParams['svg.fonttype'] = 'none'

for table in ld.TABLES.glob("*.pkl"):
    ind = table.stem
    df = pd.read_pickle(table)

    saveto = ld.RESULTS / ind
    saveto.mkdir(exist_ok=True, parents=True)

    loopsize = {}
    for category in [
        'HSV-1 induced & Z22 enriched', 'IAV induced & Z22 enriched',
        'HSV-1 induced & FLAG enriched', 'IAV induced & FLAG enriched',
    ]:
        mask = df[category]
        loopsize[category] = df.loc[mask, 'loop-size'].copy()

    loopsize = pd.DataFrame(loopsize)
    loopsize = loopsize.reset_index().melt(
        id_vars=['partition', 'coordinates'], var_name="dataset", value_name='distance'
    ).dropna()
    loopsize = loopsize[loopsize['distance'] >= 0].copy()

    maxsize = ld.c.loopsize.maxsize
    loopsize.loc[loopsize['distance'] > maxsize, 'distance'] = maxsize

    loopsize['antibody'] = loopsize['dataset'].apply(lambda x: x.split()[-2])
    loopsize['virus'] = loopsize['dataset'].apply(lambda x: x.split()[0])

    assert set(loopsize['virus']) <= {'HSV-1', 'IAV'} and set(loopsize['antibody']) == {'Z22', 'FLAG'}

    displot = sns.displot(
        loopsize, x='distance', row='virus', col='antibody',
        row_order=['HSV-1', 'IAV'], col_order=['Z22', 'FLAG'], height=2, aspect=1.5,
        kind='hist', stat='percent', facet_kws={"sharey": False}, kde=False, common_norm=False,
        binwidth=ld.c.loopsize.binwidth, binrange=(0, maxsize),
    )
    displot.fig.suptitle(ind, y=1.0)


    def annotate(data, **_kws):
        median, q90 = data['distance'].quantile([0.5, 0.9])

        ax = plt.gca()
        offset = ld.c.loopsize.binwidth / 2

        ax.axvline(median, color='black', linestyle='-')
        ax.text(
            median + offset, 1.0, f"50%\n{median / 1000:.1f} knt", rotation=0, va='top', ha='left',
            transform=ax.get_xaxis_transform(), fontsize=8
        )

        ax.axvline(q90, color='black', linestyle='--')
        ax.text(
            q90 + offset, 1.0, f"90%\n{q90 / 1000:.1f} knt", rotation=0, va='top', ha='left',
            transform=ax.get_xaxis_transform(), fontsize=8
        )

        ax.set_xlim(0, maxsize)
        xticks = list(range(0, maxsize + 1, 10_000))
        ax.set_xticks(xticks)
        ax.set_xticklabels([f"{x / 1_000:.0f} knt" for x in xticks])
        ax.set_xlabel("Distance between dsRNA arms")

        ax.set_ylim(0, 40)
        ax.set_yticks([0, 40])
        ax.set_ylabel("Proportion (%)")

        ax.text(
            1.0, 1.0, f"N={len(data)}", transform=ax.transAxes,
            verticalalignment='top', horizontalalignment='right', fontsize=12
        )


    displot.map_dataframe(annotate)

    displot.fig.savefig(saveto / "loop-size-distribution.svg", bbox_inches="tight", pad_inches=0, transparent=True)
    # displot.fig.show()
    plt.close(displot.fig)
