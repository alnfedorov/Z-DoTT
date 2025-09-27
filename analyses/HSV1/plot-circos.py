import pickle
from pathlib import Path

import numpy as np
import numpy.typing as npt
from biobit.core.loc import PerStrand
from joblib import Parallel, delayed
from matplotlib import pyplot as plt
from pycirclize import Circos
from pycirclize.track import Track

import ld
from assemblies import HSV1
from ld.enrichment import EnrichmentTracks

plt.rcParams['svg.fonttype'] = 'none'

SAVETO = ld.PLOTS / "circos"
SAVETO.mkdir(parents=True, exist_ok=True)

with open(ld.annotation.circos, 'rb') as stream:
    ANNOTATION = pickle.load(stream)


def plot_genes(track: Track, strand: str):
    tss_height = 0.035
    if strand == "+":
        vmin, vmax = 0 - tss_height, 1 + tss_height
    else:
        vmin, vmax = -1 - tss_height, 0 + tss_height

    for element, h in ("exons", tss_height), ("introns", tss_height * 0.35):
        for start, end, _, score, color in ANNOTATION[element, strand]:
            track.fill_between(
                [start, end], [score - h, score - h], [score + h, score + h],
                fc=color, vmin=vmin, vmax=vmax, lw=0
            )

    length = 1.35
    for start, end, _, score, _ in ANNOTATION["5`end", strand]:
        if strand == "+":
            r = track.r_plot_lim[0] + score * track.r_plot_size
        else:
            r = track.r_plot_lim[1] + score * track.r_plot_size
        angle = length / r
        gap = angle / track.rad_size * track.size

        if strand == "+":
            end = min(end, start + gap)
        else:
            start = max(start, end - gap)

        track.fill_between(
            [start, end], [score - tss_height, score - tss_height], [score + tss_height, score + tss_height],
            fc='black', vmin=vmin, vmax=vmax, lw=0, zorder=100
        )


def job(
        fwd: npt.NDArray[np.float32], rev: npt.NDArray[np.float32], title: str, saveto: Path
):
    circos = Circos(HSV1.segments, endspace=False, start=-270, end=90)
    circos.line(r=(45, 110), deg_lim=(90, 90), lw=4, color='black')
    sector = circos.sectors[0]

    # Outline
    track = sector.add_track(r_lim=(45, 100))
    track.axis(fc="none", lw=2, ec="black")
    track.line([0, HSV1.size], [0, 0], lw=2, color='black')

    # Ticks
    step, formatter = 25_000, lambda v: f"{v // 1000}"
    ticks = list(range(0, HSV1.size, step))
    if ticks[-1] != HSV1.size:
        ticks.append(track.size)
    labels = [formatter(v) for v in ticks]
    labels[-1] = f"{track.size / 1000:.1f} knt"
    labels[0] = "0"

    track.xticks(
        ticks, labels, tick_length=1.75, outer=True, label_orientation='horizontal', label_size=16,
        line_kws={"lw": 2}
    )

    # Terminal repeats
    for start, end in HSV1.terminal_repeats:
        track.fill_between(
            [start, end], [0, 0], [1, 1],
            fc=ld.palette.terminal, ec=None, lw=0, vmin=0, vmax=1, alpha=0.5
        )

    # Forward / Reverse genes
    track = sector.add_track((80, 100), r_pad_ratio=0.2)
    plot_genes(track, "+")

    track = sector.add_track((45, 65), r_pad_ratio=0.2)
    plot_genes(track, "-")

    # Bin the data for visualization
    # Required, because otherwise a very large SVG is produced
    cutoff = (HSV1.size // ld.BINSIZE) * ld.BINSIZE
    fwd = fwd[:cutoff].reshape(-1, ld.BINSIZE).mean(axis=1)
    rev = rev[:cutoff].reshape(-1, ld.BINSIZE).mean(axis=1)

    x = np.arange(len(fwd)) * ld.BINSIZE + ld.BINSIZE / 2

    minval = 0
    maxval = round(np.quantile(np.concatenate([fwd, rev]), 0.995), 2)

    fwd, rev = fwd.round(2), rev.round(2)
    fwd, rev = np.clip(fwd, minval, maxval), np.clip(rev, minval, maxval)

    # This is a workaround for the internal bug in the pycirclize library
    maxval *= 1.01

    # Plot the data
    track = sector.add_track((48.5, 72.5))
    mask = rev > minval
    if mask.any():
        track.bar(
            x[mask], -(rev[mask] - minval), width=ld.BINSIZE, bottom=-minval, align='center',
            vmin=-maxval, vmax=-minval, fc=ld.palette.reverse, ec=ld.palette.reverse, lw=0
        )

    track = sector.add_track((72.5, 96.5))
    mask = fwd > minval
    if mask.any():
        track.bar(
            x[mask], (fwd[mask] - minval), width=ld.BINSIZE, bottom=minval, align='center',
            vmin=minval, vmax=maxval, fc=ld.palette.forward, ec=ld.palette.reverse, lw=0
        )

    # Central line between +/- strands
    track.line([0, HSV1.size], [0, 0], lw=2, color='black')

    circos.text(f"{title}\n[{minval:.1f};{maxval:.1f}]", size=30, weight='bold')
    fig = circos.plotfig()
    fig.savefig(saveto, transparent=True, bbox_inches="tight", pad_inches=0)
    # fig.show()
    plt.close(fig)


workload = []
for pkl in ld.ENRICHMENT.glob("*.pkl"):
    with open(pkl, 'rb') as stream:
        tracks: EnrichmentTracks = pickle.load(stream)

    selected = PerStrand(forward=None, reverse=None)
    for strand in "+", "-":
        scaled = tracks.scaled_enrichment[strand][HSV1.contig]

        # Zero irrelevant values
        lowfc = tracks.enrichment[strand][HSV1.contig] < ld.MIN_ENRICHMENT
        scaled[lowfc] = 0.0
        for start, end in HSV1.terminal_repeats:
            scaled[start:end] = 0

        selected[strand] = scaled

    if "vs" in pkl.stem:
        cells, _, rip, replica = pkl.stem.split(" vs ")[0].split("+")[1].split("_")
    else:
        cells, _, rip = pkl.stem.split("_")
        replica = "merged"

    workload.append((
        selected.forward, selected.reverse,
        f"{cells}\n{rip} {replica}",
        SAVETO / pkl.with_suffix(".svg").name,
    ))

Parallel(n_jobs=-1)(delayed(job)(*args) for args in workload)
