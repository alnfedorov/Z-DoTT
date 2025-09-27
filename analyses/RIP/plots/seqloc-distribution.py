from collections import Counter

import matplotlib.pyplot as plt
import pandas as pd

import ld

plt.rcParams['svg.fonttype'] = 'none'

for table in ld.TABLES.glob("*.pkl"):
    ind = table.stem
    df = pd.read_pickle(table)

    saveto = ld.RESULTS / ind
    saveto.mkdir(exist_ok=True, parents=True)

    # Note - all sequences without inv. rev. comp. mates are classified as "Unresolved", therefore it's safe to
    # use Inv. other as a catch-all category for rare sequences.
    df['sequence[cls]'] = df['sequence[cls]'].apply(lambda x: x if x in ld.c.sequences.cls.palette else "Inv. other")
    df.loc[df['sequence[cls]'] == 'Inv. other', 'sequence[family]'] = 'Inv. other'

    for category in [
        'HSV-1 induced & Z22 enriched', 'IAV induced & Z22 enriched',
        'HSV-1 induced & FLAG enriched', 'IAV induced & FLAG enriched',
    ]:
        mask = df[category]
        N = mask.sum()
        if N == 0:
            continue

        for column, palette, order, startangle in [
            ('location', ld.c.locations.palette, ld.c.locations.order[::-1], 0),
            ('sequence[cls]', ld.c.sequences.cls.palette, ld.c.sequences.cls.order[ind][::-1], 90),
        ]:
            cnts = df.loc[mask, column].value_counts()
            cnts /= cnts.sum()

            cnts = [(x, cnts[x]) for x in order if x in cnts]
            labels, values = zip(*cnts)

            fig, ax = plt.subplots(figsize=(8, 8))
            ax.pie(
                values,
                labels=labels,
                autopct=lambda x: f"{x:.1f}%",
                colors=[palette[x] for x in labels],
                startangle=startangle, counterclock=False,
                wedgeprops=dict(width=0.5, linewidth=0.9, edgecolor='black')
            )
            ax.set_title(category + f"\nN={N}")
            fig.savefig(
                saveto / f"{category} {column}-distribution.svg",
                bbox_inches="tight", pad_inches=0, transparent=True
            )
            # plt.show()
            plt.close(fig)

        # Make a summary table for repeat families
        cnts = Counter(df.loc[mask, ['sequence[cls]', 'sequence[family]']].itertuples(index=False, name=None))
        cnts = pd.Series(cnts).to_frame('N').reset_index(names=['Class', 'Family'])
        cnts = cnts.sort_values(['Class', 'N'], ascending=False).set_index('Class')
        cnts.to_csv(saveto / f"{category} sequences-summary.csv")
