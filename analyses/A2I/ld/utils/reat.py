import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

REAT = "reat"


def _reat(command: list[str], buffer: Path) -> pd.DataFrame:
    command.extend(["-o", buffer])
    subprocess.run(command, check=True)
    assert buffer.exists()

    sites = pd.read_csv(buffer, dtype={
        'contig': str, 'refnuc': str, 'prednuc': str, 'trstrand': str,
        'pos': np.uint32, "A": np.uint32, "C": np.uint32, "G": np.uint32, "T": np.uint32
    })
    buffer.unlink()
    return sites


def run(
        fasta: Path, inputs: list[Path], stranding: str, paired: bool, saveto: Path, exclude: Path,
        min_coverage: int = 6, min_edits: int = 1, freqthr: tuple[float, float] = (0.01, 0.50), threads: int = 64
):
    assert saveto.suffix == '.csv'
    if saveto.with_suffix(".csv.gz").exists():
        return

    command = [
        REAT, "site", "-i", *inputs, "--exclude", exclude,
        "-r", fasta, "--stranding", stranding, "-t", str(threads),
        "--mapq", "254", "--ex-flags", "2828", "--phread", "25",
        "--out-min-cov", str(min_coverage), "--out-min-freq", '0', "--out-min-mismatches", str(min_edits)
    ]
    if paired:
        command.extend(["--in-flags", "3"])

    # Gather candidates & keep only A-to-I edits
    esites = _reat(command, saveto)

    fwd = (esites['trstrand'] == "+") & \
          (esites['A'] > 0) & \
          (esites['refnuc'] == 'A') & \
          ((esites['A'] + esites['G']) > (esites['T'] + esites['C']))
    rev = (esites['trstrand'] == "-") & \
          (esites['T'] > 0) & \
          (esites['refnuc'] == 'T') & \
          ((esites['T'] + esites['C']) > (esites['A'] + esites['G']))
    esites = esites[fwd | rev].copy()

    esites['refnuc'] = np.where(esites['trstrand'] == "+", esites['A'], esites['T'])
    esites['misnuc'] = np.where(esites['trstrand'] == "+", esites['G'], esites['C'])

    freq = esites['misnuc'] / (esites['refnuc'] + esites['misnuc'])
    esites = esites[(freq >= freqthr[0]) & (freq <= freqthr[1])]

    esites = esites[(esites['refnuc'] + esites['misnuc']) >= min_coverage]
    esites = esites[esites['misnuc'] >= min_edits]

    esites = esites[['contig', 'trstrand', 'pos', 'refnuc', 'misnuc']]

    saveto = saveto.with_suffix(".csv.gz")
    esites.to_csv(saveto, index=False)
