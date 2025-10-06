from pathlib import Path
from subprocess import check_call
from tempfile import NamedTemporaryFile
from typing import Literal, Any

import pandas as pd
from attrs import define
from biobit.toolkit.libnorm import MedianOfRatiosNormalization
from joblib import Parallel, delayed

_DE_ANALYSIS_WRAPPER = Path(__file__).parent / "deanalysis.R"

_ARGUMENTS_SEPARATOR = "_-$-_"


@define(slots=True, frozen=True, eq=True, order=True, repr=True, hash=True)
class Comparison:
    # Unique identifier
    ind: str
    # Sample IDs must match column names in the counts table
    samples: tuple[Any, ...]
    # Design formula
    design: str
    # Alternative hypothesis
    alternative: Literal["greaterAbs", "lessAbs", "greater", "less"]
    # Attribute to test
    attribute: str
    # Target and reference levels
    target: str
    reference: str
    # Thresholds
    log2fc_thr: float
    padj_thr: float
    min_counts: int
    min_replicates: int

    def __attrs_post_init__(self):
        assert self.target != self.reference, "Target and reference levels must be different."

    def normalize(self, smplmap: dict[Any, str], attrmap: dict[str, str], lvlmap: dict[str, str]) -> 'Comparison':
        design = self.design
        for old, new in attrmap.items():
            design = design.replace(old, new)
        samples = tuple(smplmap[smpl] for smpl in self.samples)

        return Comparison(
            self.ind, samples, design, self.alternative,
            attrmap[self.attribute], lvlmap[self.target], lvlmap[self.reference],
            self.log2fc_thr, self.padj_thr, self.min_counts, self.min_replicates
        )


@define(slots=True, frozen=True, eq=True, order=True, repr=True, hash=True)
class DESeq2:
    counts: pd.DataFrame
    samples: pd.DataFrame
    normalization: MedianOfRatiosNormalization

    comparisons: tuple[Comparison, ...]

    def __attrs_post_init__(self):
        samples = set(self.samples.index)
        if samples != set(self.counts.columns):
            raise ValueError("Samples index and counts table columns must be identical.")

        inds = set()
        for cmp in self.comparisons:
            assert cmp.ind not in inds, f"Duplicate comparison identifier {cmp.ind}"
            inds.add(cmp.ind)

            assert all(smpl in samples for smpl in cmp.samples), f"Unknown samples in comparison {cmp}"
            assert cmp.attribute in self.samples.columns, f"Unknown attribute {cmp.attribute} in comparison {cmp}"

            levels = set(self.samples[cmp.attribute].unique())
            assert cmp.target in levels, f"Unknown target level {cmp.target} in comparison {cmp}"
            assert cmp.reference in levels, f"Unknown reference level {cmp.reference} in comparison {cmp}"

            assert 0 <= cmp.log2fc_thr, f"Log2FC threshold must be non-negative, got {cmp.log2fc_thr}"
            assert 0 <= cmp.padj_thr <= 1, f"Adjusted p-value threshold must be in the range [0, 1], got {cmp.padj_thr}"

    def _normalize_string(self, string: str) -> str:
        allsymbols = set(string)
        for symbol in allsymbols:
            if symbol.isalnum() or symbol == "_":
                continue
            string = string.replace(symbol, "_")

        # Remove multiple underscores
        string = "_".join([x for x in string.split("_") if x])
        return string

    def run(self, backend: Parallel) -> dict[str, pd.DataFrame]:
        # Map all samples to simple ids like: sample_1, sample_2, etc.
        smplmap = {smpl: self._normalize_string(str(smpl)) + f"_{ind}" for ind, smpl in enumerate(self.samples.index)}
        samples = self.samples.rename(index=smplmap)
        counts = self.counts.rename(columns=smplmap)
        normalization = self.normalization.rename(smplmap)

        # Map all attributes and levels to simple ids like: attribute_1, attribute_2, etc.
        attrmap = {attr: self._normalize_string(str(attr)) + f"_{ind}" for ind, attr in enumerate(self.samples.columns)}
        samples = samples.rename(columns=attrmap)

        all_lvls = set(x for x in samples.values.flatten())
        lvlmap = {lvl: self._normalize_string(str(lvl)) + f"_{ind}" for ind, lvl in enumerate(all_lvls)}

        for attr in samples.columns:
            samples[attr] = samples[attr].map(lvlmap)

        # Remap counts index to plain ids
        featmap = {feat: self._normalize_string(str(feat)) + f"_{ind}" for ind, feat in enumerate(counts.index)}
        counts = counts.rename(index=featmap)

        # Launch the DESeq2 analysis for each comparison
        deseq2 = backend(
            delayed(_run)(counts, samples, normalization, cmp.normalize(smplmap, attrmap, lvlmap))
            for cmp in self.comparisons
        )

        # Postprocess the results
        featrevmap = {v: k for k, v in featmap.items()}
        result = {cmpind: df.rename(index=featrevmap) for cmpind, df in deseq2}

        return result


def _run(counts: pd.DataFrame, samples: pd.DataFrame, normalization: MedianOfRatiosNormalization, cmp: Comparison):
    # Subsample the counts table to include only the samples of interest and normalize it
    counts = counts[list(cmp.samples)].copy()
    _, _, scaling = normalization.scaling_factors({x: [x] for x in cmp.samples})
    for key, val in scaling.to_dict().items():
        counts[key] = (counts[key] / val).round(0).astype(int)

    # Subsample to include only regions with enough counts
    counts = counts.loc[(counts >= cmp.min_counts).sum(axis=1) >= cmp.min_replicates].copy()

    # Subsample the samples table to include only the samples of interest
    samples = samples.loc[list(cmp.samples)]

    assert cmp.attribute in samples.columns
    assert cmp.target in samples[cmp.attribute].unique()
    assert cmp.reference in samples[cmp.attribute].unique()

    # Setup and run the DESeq2
    with (
        NamedTemporaryFile() as smpltable,
        NamedTemporaryFile() as cntstable,
        NamedTemporaryFile() as saveto
    ):
        samples.to_csv(smpltable, index=True)
        counts.to_csv(cntstable, index=True)

        parameters = [
            "Rscript", _DE_ANALYSIS_WRAPPER, smpltable.name, cntstable.name, cmp.design,
            cmp.attribute, cmp.target, cmp.reference, cmp.alternative,
            str(cmp.log2fc_thr), str(cmp.padj_thr), saveto.name,
        ]
        print(parameters, "\n")
        check_call(parameters)
        results = pd.read_csv(saveto, compression="gzip", index_col=0)
    return cmp.ind, results
