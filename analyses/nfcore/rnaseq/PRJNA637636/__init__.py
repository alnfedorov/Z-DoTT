# AUTO-GENERATED FILE - DO NOT EDIT
from pathlib import Path

from biobit.toolkit import nfcore, seqproj

ROOT = Path(__file__).parent


def PRJNA637636() -> seqproj.Project:
    from data import PRJNA637636
    project = PRJNA637636()
    
    return nfcore.rnaseq.parse.into_seqproj(
        project, ROOT / 'results',
        seqexp2descriptor=lambda exp: nfcore.rnaseq.descriptor.from_seqexp(exp, title_builder="title")
    )
