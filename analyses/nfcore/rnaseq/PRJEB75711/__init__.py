# AUTO-GENERATED FILE - DO NOT EDIT
from pathlib import Path

from biobit.toolkit import nfcore, seqproj

ROOT = Path(__file__).parent


def PRJEB75711() -> seqproj.Project:
    from data import PRJEB75711
    project = PRJEB75711()
    
    return nfcore.rnaseq.parse.into_seqproj(
        project, ROOT / 'results',
        seqexp2descriptor=lambda exp: nfcore.rnaseq.descriptor.from_seqexp(exp, title_builder="title")
    )
