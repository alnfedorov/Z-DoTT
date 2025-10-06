from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.transforms import blended_transform_factory

plt.rcParams['svg.fonttype'] = 'none'


class VolcanoPreset:
    def __init__(self, data, log10pv: str, padj: str, log2fc: str):
        self.payload = pd.DataFrame({
            "log2fc": data[log2fc], "-log10pv": data[log10pv], "padj": data[padj],
        })
        self.payload['color'] = None
        self.payload['alpha'] = None
        self.payload['lw'] = None
        self.payload['annotation'] = None

        self.data = data

        self.padj = None
        self.log2fc = None
        self.legend = None
        self.xlim = (None, None)
        self.ylim = (None, None)

    def color(self, column: str, mapping: dict[Any, str]):
        self.payload['color'] = self.data[column].map(mapping)
        return self

    def alpha(self, column: str, mapping: dict[Any, float]):
        self.payload['alpha'] = self.data[column].map(mapping)
        return self

    def edge_width(self, column: str, mapping: dict[Any, float]):
        self.payload['lw'] = self.data[column].map(mapping)
        return self

    def size(self, column: str, mapping: dict[Any, float]):
        self.payload['s'] = self.data[column].map(mapping)
        return self

    def annotate(self, column: str, mapping: dict[Any, str]):
        self.payload['annotation'] = self.data[column].map(mapping)
        return self

    def marker(self, column: str, mapping: dict[Any, str]):
        self.payload['marker'] = self.data[column].map(mapping)
        return self

    def threshold(self, padj: float | None = None, log2fc: float | None = None):
        self.padj = padj
        self.log2fc = log2fc
        return self

    def xlimit(self, xmin: float | None = None, xmax: float | None = None):
        self.xlim = (xmin, xmax)
        return self

    def ylimit(self, ymin: float | None = None, ymax: float | None = None):
        self.ylim = (ymin, ymax)
        return self

    def set(
            self,
            size: float | None = None,
            alpha: float | None = None,
            annotate: str | None = None,
            line_width: float | None = None,
    ):
        if size is not None:
            self.payload['s'] = size
        if alpha is not None:
            self.payload['alpha'] = alpha
        if annotate is not None:
            self.payload['annotation'] = self.data[annotate]
        if line_width is not None:
            self.payload['lw'] = line_width
        return self

    def plot(self, ax: plt.Axes, **kwargs) -> plt.Axes:
        # Make markers to annotate dots outside the limits
        if 'marker' not in self.payload:
            self.payload['marker'] = 'o'

        if self.ylim[0] is not None:
            mask = self.payload['-log10pv'] < self.ylim[0]
            self.payload.loc[mask, 'marker'] = 'v'
            self.payload.loc[mask, '-log10pv'] = self.ylim[0]

        if self.ylim[1] is not None:
            mask = self.payload['-log10pv'] > self.ylim[1]
            self.payload.loc[mask, 'marker'] = '^'
            self.payload.loc[mask, '-log10pv'] = self.ylim[1]

        if self.xlim[0] is not None:
            mask = self.payload['log2fc'] < self.xlim[0]
            self.payload.loc[mask, 'marker'] = '<'
            self.payload.loc[mask, 'log2fc'] = self.xlim[0]

        if self.xlim[1] is not None:
            mask = self.payload['log2fc'] > self.xlim[1]
            self.payload.loc[mask, 'marker'] = '>'
            self.payload.loc[mask, 'log2fc'] = self.xlim[1]

        # Use default values for missing parameters
        keys = ['color', 'alpha', 'lw', 'marker', 's']
        for k in keys:
            if k not in self.payload.columns:
                self.payload[k] = None

        # Plot the data
        for config, data in self.payload.groupby(keys, dropna=False):
            params = {'edgecolor': 'black'}
            for k, v in zip(keys, config):
                if not pd.isna(v):
                    params[k] = v
            ax.scatter(data['log2fc'], data['-log10pv'], **(kwargs | params))

        # Set the limits
        limits = ax.get_xlim()
        if self.xlim[0] is not None:
            limits = self.xlim[0], limits[1]
        if self.xlim[1] is not None:
            limits = limits[0], self.xlim[1]
        ax.set_xlim(limits)

        limits = ax.get_ylim()
        if self.ylim[0] is not None:
            limits = self.ylim[0], limits[1]
        if self.ylim[1] is not None:
            limits = limits[0], self.ylim[1]
        ax.set_ylim(limits)

        ax.set(xlabel="$log_{2}$(fold-change)", ylabel="$-log_{10}$(p-value)")

        # Plot thresholds
        if self.log2fc is not None:
            ax.axvline(self.log2fc, ls='--', color='black')
            ax.axvline(-self.log2fc, ls='--', color='black')
            # transform = blended_transform_factory(ax.transData, ax.transAxes)
            # ax.text(0, 1.0, f"$\\bf{{FC~≥~{2 ** self.log2fc:.2f}}}$", ha='center', va='bottom', transform=transform)

        if self.padj is not None:
            significant = self.payload[self.payload['padj'] <= self.padj]
            if len(significant) > 0:
                threshold = significant['-log10pv'].min()
                ax.axhline(threshold, ls='--', color='black')
                transform = blended_transform_factory(ax.transAxes, ax.transData)
                ax.text(1.0, threshold, f'$\\bf{{FDR~≤~{self.padj}}}$', ha='right', va='bottom', transform=transform)

        # Annotate hits
        if 'annotation' in self.payload:
            hits = self.payload[~self.payload['annotation'].isna()]
            for _, h in hits.iterrows():
                x, y, name = h['log2fc'], h['-log10pv'], h['annotation']
                if h['padj'] <= self.padj and abs(h['log2fc']) >= self.log2fc:
                    ax.text(x, y, name, fontsize=8, fontdict={'fontweight': 'bold'})
                else:
                    ax.text(x, y, name, fontsize=6)
        return ax
