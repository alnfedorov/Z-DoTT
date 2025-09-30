# AUTO-GENERATED FILE - DO NOT EDIT
from pathlib import Path

from biobit.toolkit import nfcore, seqproj

ROOT = Path(__file__).parent


def PRJNA382632() -> seqproj.Project:
    from data import PRJNA382632
    return nfcore.rnaseq.parse.into_seqproj(
        PRJNA382632(), ROOT / 'results',
        seqexp2descriptor=lambda exp: nfcore.rnaseq.descriptor.from_seqexp(exp, title_builder="title")
    )
