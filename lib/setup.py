import subprocess
from pathlib import Path
from subprocess import check_call
from typing import Literal


def download(url: str, saveto: Path):
    """Download a file from a URL to a specified path"""
    saveto.parent.mkdir(parents=True, exist_ok=True)
    check_call(['wget', url, '-O', saveto])


def ungzip_then_bgzip(gzfile: Path):
    """Ungzip a file and then bgzip it inplace"""
    tmp = gzfile.with_suffix('.tmp')
    gzfile.rename(tmp)
    check_call(f"gunzip {tmp} --stdout | bgzip --threads $(nproc) -l 6 -o {gzfile} /dev/stdin", shell=True)
    tmp.unlink()


def index_fasta(fasta: Path):
    """Index a FASTA file using samtools"""
    check_call(['samtools', 'faidx', fasta])


def run_in_pixi(environment: Literal["nextflow", "liftoff"], cwd: Path, cmd: str):
    subprocess.check_call(f'pixi run -e {environment} "{cmd}"', shell=True, cwd=cwd)
