# AUTO-GENERATED FILE - DO NOT EDIT
from pathlib import Path

from biobit.toolkit import seqproj

ROOT = Path(__file__).parent


def PRJNA637636() -> seqproj.Project:
    return seqproj.adapter.yaml.load(ROOT / "seq-project.yaml")
