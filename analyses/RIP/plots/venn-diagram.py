import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn2

import ld

plt.rcParams['svg.fonttype'] = 'none'

for table in ld.TABLES.glob("*.pkl"):
    ind = table.stem
    df = pd.read_pickle(table)

    saveto = ld.RESULTS / ind
    saveto.mkdir(exist_ok=True, parents=True)

    for title, (x, y) in {
        "HSV-1 vs IAV [Z22]": [
            "HSV-1 induced & Z22 enriched",
            "IAV induced & Z22 enriched",
        ],
        "HSV-1 vs IAV [FLAG]": [
            "HSV-1 induced & FLAG enriched",
            "IAV induced & FLAG enriched",
        ],
        "HSV-1 [Z22-vs-FLAG]": [
            "HSV-1 induced & Z22 enriched",
            "HSV-1 induced & FLAG enriched",
        ],
        "IAV [Z22-vs-FLAG]": [
            "IAV induced & Z22 enriched",
            "IAV induced & FLAG enriched",
        ]
    }.items():
        _x, _y = set(df[df[x]].index), set(df[df[y]].index)

        plt.figure(figsize=(8, 8))
        venn2((_x, _y), (x, y))
        plt.title(f"{ind}\n{title}", fontsize=16)
        plt.tight_layout()
        plt.savefig(saveto / f"{title} overlap.svg", bbox_inches="tight", pad_inches=0, transparent=True)
        # plt.show()
        plt.close()

    # Overlap between viruses by using the union of detected Z-RNAs
    x = set(df[df["HSV-1 induced & Z22 enriched"]].index) | set(df[df["HSV-1 induced & FLAG enriched"]].index)
    y = set(df[df["IAV induced & Z22 enriched"]].index) | set(df[df["IAV induced & FLAG enriched"]].index)

    plt.figure(figsize=(8, 8))
    venn2((x, y), ("HSV-1 induced (Z22 or FLAG)", "IAV induced (Z22 or FLAG)"))
    plt.title(f"{ind}\nHSV-1 vs IAV [Z22 or FLAG]", fontsize=16)
    plt.tight_layout()
    plt.savefig(saveto / "HSV-1 vs IAV [Z22 or FLAG] overlap.svg", bbox_inches="tight", pad_inches=0, transparent=True)
    # plt.show()
    plt.close()
