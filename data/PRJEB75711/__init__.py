# AUTO-GENERATED FILE - DO NOT EDIT
from pathlib import Path

from biobit.toolkit import seqproj

ROOT = Path(__file__).parent


def PRJEB75711() -> seqproj.Project:
    return seqproj.adapter.yaml.load(ROOT / "seq-project.yaml")
