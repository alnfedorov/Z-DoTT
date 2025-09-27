import pandas as pd
from attrs import define

from stories.aberrantome import Config


@define
class Workload:
    comparison: Config
    aberrantome: dict[str, pd.DataFrame]


palette = {
    "Not significant": "#B0B0B0",
    "Significant up": "#FFB74D",
    "Strong up": "#D32F2F",
    "Significant down": "#8BC34A",
    "Strong down": "#388E3C"
}


class threshold:
    # Thresholds for significance (dexseq-based stat test)
    padj = 0.01
    log2fc = 0.585

    # Minimum delta between aberrant transcription rates to consider the change "strong"
    delta = 0.1

    trlimit = 0.5


class annotate:
    by_assembly = {
        "CHM13v2": [
            "C6orf62", "RUVBL2", "RABGGTB", "ROCK1", "NCL", "KBTBD2", "MSH4", "CGB3", "CGB2", "ZNF132", "ZNF584",
        ],
        "GRCm39": [
            "Btbd3", "Ptbp1", "Gtpbp4", "H2ac18", "H2ac19", "Nabp1", "Hmga1", "Rock1", "Vbp1", "Gm14419",
            "Gnpda1", "Anapc4", "Pcna", "Haus2", "Tlcd1", "2310057M21Rik", "Fam24a", "Fam24b"
        ]
    }
