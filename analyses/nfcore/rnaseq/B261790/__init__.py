# AUTO-GENERATED FILE - DO NOT EDIT
from pathlib import Path

from biobit.toolkit import nfcore, seqproj

ROOT = Path(__file__).parent


def B261790() -> seqproj.Project:
    from data import B261790
    return nfcore.rnaseq.parse.into_seqproj(
        B261790(), ROOT / 'results',
        seqexp2descriptor=lambda exp: nfcore.rnaseq.descriptor.from_seqexp(exp, title_builder="title")
    )
