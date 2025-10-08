import subprocess
from pathlib import Path
from subprocess import check_call
from typing import Literal


def download(url: str, saveto: Path):
    """Download a file from a URL to a specified path"""
    saveto.parent.mkdir(parents=True, exist_ok=True)
    check_call(['wget', url, '-O', saveto])


def rebgzip(gzfile: Path):
    """Ungzip a file and then bgzip it inplace"""
    tmp = gzfile.with_suffix('.tmp')
    gzfile.rename(tmp)
    check_call(f"gunzip {tmp} --stdout | bgzip --threads $(nproc) -l 6 -o {gzfile} /dev/stdin", shell=True)
    tmp.unlink()


def bgzip(file: Path, saveto: Path):
    """Bgzip a file inplace or to a specified location"""
    saveto.parent.mkdir(parents=True, exist_ok=True)
    check_call(f"bgzip --threads $(nproc) -l 6 -c {file} > {saveto}", shell=True)


def gzip(file: Path, saveto: Path):
    """Gzip a file to a specified location"""
    saveto.parent.mkdir(parents=True, exist_ok=True)
    check_call(f"gzip -c {file} > {saveto}", shell=True)


def faidx(fasta: Path):
    """Index a FASTA file using samtools"""
    check_call(['samtools', 'faidx', fasta])


def autocat(files: list[Path], outfile: Path):
    """Concatenate multiple files into one using cat or zcat based on file extension"""
    outfile.parent.mkdir(parents=True, exist_ok=True)
    with open(outfile, 'w') as out:
        for file in files:
            cmd = 'zcat' if file.suffix in {'.gz', '.bgz'} else 'cat'
            subprocess.check_call(f"{cmd} {file}", shell=True, stdout=out)


def run_in_pixi(environment: Literal["nextflow", "liftoff"], cwd: Path, cmd: str):
    subprocess.check_call(f'pixi run -e {environment} "{cmd}"', shell=True, cwd=cwd)
