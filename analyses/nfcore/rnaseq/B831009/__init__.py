# AUTO-GENERATED FILE - DO NOT EDIT
from pathlib import Path

from biobit.toolkit import nfcore, seqproj

ROOT = Path(__file__).parent


def B831009() -> seqproj.Project:
    from data import B831009
    return nfcore.rnaseq.parse.into_seqproj(
        B831009(), ROOT / 'results',
        seqexp2descriptor=lambda exp: nfcore.rnaseq.descriptor.from_seqexp(exp, title_builder="title")
    )
