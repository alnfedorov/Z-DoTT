import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from biobit.core.loc import PerStrand
from joblib import Parallel, delayed

import ld
from assemblies import HSV1
from ld.enrichment import EnrichmentTracks

plt.rcParams['svg.fonttype'] = 'none'

SAVETO = ld.PLOTS / "enrichment-ratio"
SAVETO.mkdir(parents=True, exist_ok=True)


def job(path: Path):
    with open(path, "rb") as stream:
        data: EnrichmentTracks = pickle.load(stream)

    mask = PerStrand(
        forward=data.enrichment.forward[HSV1.contig] < ld.MIN_ENRICHMENT,
        reverse=data.enrichment.reverse[HSV1.contig] < ld.MIN_ENRICHMENT
    )
    for postfix, values in ("enrichment", data.enrichment), ("scaled_enrichment", data.scaled_enrichment):
        fwd, rev = values.forward[HSV1.contig], values.reverse[HSV1.contig]
        fwd[mask.forward], rev[mask.reverse] = 0, 0

        # Calculate the difference
        diff = (fwd - rev) / (fwd + rev + 1e-12)
        for start, end in HSV1.terminal_repeats:
            diff[start:end] = 0
        # diff[:] = uniform_filter1d(diff, size=ld.BINSIZE * 3)  # Smooth the curve for better visualization

        # Make the plot
        fig, axes = plt.subplots(2, 1, sharex=True, sharey=False, figsize=(12, 6))

        cutoff = (HSV1.size // ld.BINSIZE) * ld.BINSIZE
        X = np.arange(0, HSV1.size)[:cutoff].reshape(-1, ld.BINSIZE).mean(axis=-1) + ld.BINSIZE / 2

        # Raw enrichment values
        fwd = fwd[:cutoff].reshape(-1, ld.BINSIZE).mean(axis=-1)
        axes[0].plot(X, fwd, lw=0.5, color=ld.palette.forward)

        rev = rev[:cutoff].reshape(-1, ld.BINSIZE).mean(axis=-1)
        axes[0].plot(X, -rev, lw=0.5, color=ld.palette.reverse)

        lims = axes[0].get_ylim()
        lims = max(abs(lims[0]), abs(lims[1]))
        axes[0].set_ylim(-lims, lims)

        # Ratio
        _diff = diff[:cutoff].reshape(-1, ld.BINSIZE).mean(axis=-1)
        axes[1].plot(X, _diff, color='black', lw=0.5)

        axes[1].axhline(ld.MAX_ENRICHMENT_DIFF, lw=1, color='black', ls='--')
        axes[1].axhline(-ld.MAX_ENRICHMENT_DIFF, lw=1, color='black', ls='--')

        axes[1].set_ylim(-1, 1)

        rois = (diff <= ld.MAX_ENRICHMENT_DIFF) & (diff >= -ld.MAX_ENRICHMENT_DIFF)
        rois = np.concatenate(([0], rois, [0]))
        rois = np.abs(np.diff(rois))
        rois = np.where(rois == 1)[0].reshape(-1, 2)

        for start, end in rois[1:-1]:
            if end - start > 150:
                axes[1].axvspan(start, end, alpha=0.5, color='green')

        for ax in axes:
            for start, end in HSV1.terminal_repeats:
                ax.axvspan(start, end, alpha=0.5, color='gray')

            ax.set_xlim(0, HSV1.size)
            ax.axhline(0, lw=1, color='black')

        fig.suptitle(f"{path.name} - {postfix}")
        fig.tight_layout()
        saveto = SAVETO / path.with_suffix(f".{postfix}.svg").name
        fig.savefig(saveto, bbox_inches="tight", pad_inches=0, transparent=True)
        # fig.show()
        plt.close(fig)


Parallel(n_jobs=-1)(delayed(job)(path) for path in ld.ENRICHMENT.glob("*.pkl"))
