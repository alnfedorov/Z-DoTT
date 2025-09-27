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

    records = []
    for category in [
        'HSV-1 induced & Z22 enriched', 'IAV induced & Z22 enriched',
        'HSV-1 induced & FLAG enriched', 'IAV induced & FLAG enriched',
    ]:
        subdf = df[df[category]]
        records.append({
            "category": category,
            "N": len(subdf),
            "A-I edited (%)": subdf['A-I edited'].mean() * 100,
            "A-I edited": subdf['A-I edited'].sum(),
        })
    df = pd.DataFrame(records)

    df['virus'] = df['category'].apply(lambda x: x.split()[0])
    df['antibody'] = df['category'].apply(lambda x: x.split()[3])

    order = ['HSV-1', 'IAV']
    hue_order = ['Z22', 'FLAG']

    fig, ax = plt.subplots(figsize=(4, 8))
    sns.barplot(
        df, x='virus', y='A-I edited (%)', hue='antibody', ax=ax, lw=1.0, edgecolor='black',
        order=order, hue_order=hue_order
    )

    pos = -0.2
    for virus in order:
        for antibody in hue_order:
            record = df[(df['virus'] == virus) & (df['antibody'] == antibody)].iloc[0]
            value = record['A-I edited']
            y = record['A-I edited (%)']
            ax.text(pos, y, f"N={value:,}", ha='center', va='bottom')
            pos += 0.5
        pos + 0.2

    ax.set_title(ind)
    ax.set_ylim(0, 100)
    ax.spines[['top', 'right']].set_visible(False)

    fig.savefig(saveto / "A-I editing summary.svg", bbox_inches="tight", pad_inches=0, transparent=True)
