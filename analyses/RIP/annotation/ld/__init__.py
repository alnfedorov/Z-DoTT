import pickle
from pathlib import Path

import pandas as pd
from biobit.core.loc import Orientation, Interval

from .config import Config, StatTest

ROOT = Path(__file__).parent

RESULTS = ROOT / "results"


class cache:
    class signal:
        pkl = RESULTS / "signal.pkl"

        @staticmethod
        def load() -> dict[str, dict[tuple[str, Orientation], dict[tuple[str, str], dict[Interval, float]]]]:
            with open(cache.signal.pkl, 'rb') as stream:
                return pickle.load(stream)

    class counts:
        pkl = RESULTS / "counts.pkl"

        @staticmethod
        def load() -> dict[str, pd.DataFrame]:
            with open(cache.counts.pkl, 'rb') as stream:
                return pickle.load(stream)

    class regions:
        root = ROOT / "regions"

        @staticmethod
        def load(assembly: str) -> pd.DataFrame:
            return pd.read_csv(cache.regions.root / f"{assembly}.tsv", sep='\t', dtype={
                "seqid": str, 'Ensembl ID': str
            })

        @staticmethod
        def save(assembly: str, df: pd.DataFrame):
            df.to_csv(cache.regions.root / f"{assembly}.tsv", sep='\t', index=False)
