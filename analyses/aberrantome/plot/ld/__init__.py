from pathlib import Path

from . import config as c
from .resolve import resolve

ROOT = Path(__file__).parent
RESULTS = ROOT / "results"
WORKLOAD = RESULTS / "payload.pkl"
