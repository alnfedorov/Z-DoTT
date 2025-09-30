# AUTO-GENERATED FILE - DO NOT EDIT
from pathlib import Path

from biobit.toolkit import nfcore, seqproj

ROOT = Path(__file__).parent


def B256178() -> seqproj.Project:
    from data import B256178
    project = B256178()
    
    return nfcore.rnaseq.parse.into_seqproj(
        project, ROOT / 'results',
        seqexp2descriptor=lambda exp: nfcore.rnaseq.descriptor.from_seqexp(exp, title_builder="title")
    )
