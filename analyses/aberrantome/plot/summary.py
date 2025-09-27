import matplotlib.pyplot as plt
import pandas as pd

import ld

plt.rcParams['svg.fonttype'] = 'none'

for ind, wd in pd.read_pickle(ld.WORKLOAD).items():  # type: str, ld.config.Workload
    data = []
    asessed_genes = {}
    for name, df in wd.aberrantome.items():
        record = df['category'].value_counts().to_dict()
        total = sum(record.values())
        asessed_genes[name] = total
        assert total > 0

        record = {k: v / total * 100 for k, v in record.items()}
        record['Test'] = name
        data.append(record)

    data = pd.DataFrame(data).fillna(0)
    data = data.melt(id_vars=['Test'], var_name='Category', value_name='Percentage')
    data = data[data['Category'] != 'Not significant'].copy()

    assert data['Test'].notnull().all(), data

    fig, ax = plt.subplots(figsize=(6, 7))

    order = ['Intronic', 'Divergent', 'Upstream', 'Downstream']
    hue_order = ['Strong down', 'Significant down', 'Significant up', 'Strong up', ]
    assert set(data['Category']) <= set(hue_order), data['Category'].unique()
    assert set(data['Test']) <= set(order), data['Test'].unique()

    # Stacked bar plot
    counts = data.set_index(['Test', 'Category'])['Percentage'].to_dict()
    for i, test in enumerate(order):
        bottom = 0
        for j, category in enumerate(hue_order):
            if (test, category) not in counts or counts[test, category] == 0:
                continue

            value = counts[test, category]
            ax.bar(i, value, bottom=bottom, color=ld.c.palette[category], edgecolor='black', align='center')
            ax.text(i, bottom + value / 2, f"{value:.1f}%", ha='center', va='center', fontsize=12, color='black')
            bottom += value

    # Axes
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlim(-0.5, len(order) - 0.5)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)

    ax.set_xlabel("Aberrant transcription type", fontsize=12)
    ax.set_ylim(0, 100)

    # Total number of genes in each comparison
    for x, lbl in enumerate(order):
        ax.text(
            x, ax.get_ylim()[1], f"N={asessed_genes[lbl]:,}", ha='center', va='bottom', fontweight='bold', fontsize=8
        )

    # Legend
    legend = []
    for key in hue_order[::-1]:
        legend.append(plt.Rectangle((0, 0), 1, 1, color=ld.c.palette[key], label=key))
    ax.legend(handles=legend, loc='upper left', fontsize=8, title=None, title_fontsize=8, frameon=False)

    ax.set_title(ind, fontsize=16, fontweight='bold', y=1.025)

    # Save the results
    # fig.show()

    saveto = ld.RESULTS / ind / "summary.svg"
    saveto.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(saveto, bbox_inches="tight", pad_inches=0, transparent=True, dpi=400)

    plt.close(fig)
