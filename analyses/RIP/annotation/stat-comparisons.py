import os
import pickle
import tempfile
from collections import defaultdict
from pathlib import Path
from subprocess import check_call

import pandas as pd
from joblib import Parallel, delayed

import ld
from stories import normalization
from stories.RIP import annotation

# Load cached counts
COUNTS = ld.cache.counts.load()

# Load cached table for scaling
SCALING = normalization.median_of_ratios()

DESEQ2 = Path(__file__).parent / "deseq2.R"
assert DESEQ2.is_file()


def job(config: str, host: str, test: annotation.StatTest):
    samples = sorted(test.design.index)

    # Subsample the counts table
    counts = COUNTS[config][samples]
    assert counts.isna().sum().sum() == 0

    # Calculate the scaling for all samples
    _, _, scaling = SCALING[host].scaling_factors({x: [x] for x in samples})

    # Remap the names of samples to please DESeq2
    mapping = {k: f"sample_{ind}" for ind, k in enumerate(samples)}
    order = [mapping[x] for x in samples]

    counts = counts.rename(columns=mapping)
    counts.index.name = 'index'

    scaling.index = [mapping[x] for x in scaling.index]
    scaling = scaling.to_frame('sizeFactor').T[order]

    design = test.design.copy()
    design.index = [mapping[x] for x in design.index]
    design = design.T[order].T

    # Order must be the same
    assert (design.index == counts.columns).all()
    assert (scaling.columns == counts.columns).all()

    # Save tables
    paths = []
    for table in counts, scaling, design:
        fd, path = tempfile.mkstemp()
        os.close(fd)
        table.to_csv(path, sep='\t', index=True)
        paths.append(path)
    cntfile, scfile, desfile = paths

    # Run deseq2
    tmpdir = tempfile.mkdtemp()

    baselines = "$".join(f"{k}%{v}" for k, v in test.baselines.items())
    comparisons = "$".join(f"{k}%{v}" for k, v in test.comparisons.items())

    command = [
        "Rscript", DESEQ2, desfile, cntfile, scfile,
        test.formula, baselines, comparisons, test.alternative, str(test.log2fc), str(test.alpha), tmpdir
    ]
    print(command)
    check_call(command)

    tables = {}
    for file in Path(tmpdir).iterdir():
        df = pd.read_csv(file, index_col=0)
        df['padj'] = df['padj'].fillna(1.0)
        df['pvalue'] = df['pvalue'].fillna(1.0)

        # # R-based FDR correction is buggy (not sure why)
        # from statsmodels.stats.multitest import multipletests
        # df['padj'] = multipletests(df['pvalue'], alpha=0.01, method='fdr_bh')[1]

        key = file.name.split(".")[0]
        tables[key] = df
        file.unlink()
    os.rmdir(tmpdir)

    assert set(tables) == set(test.comparisons), (list(tables.keys()), list(test.comparisons.keys()))
    return config, tables


# Load all the configs
confs = annotation.Config.load()

# Run tests in parallel
results = Parallel(n_jobs=-1, backend='multiprocessing', verbose=100)(
    delayed(job)(config.ind, config.comparisons[0].assembly, test)
    for config in confs for test in config.tests
)

pooled = defaultdict(dict)
for ind, tables in results:
    pooled[ind].update(tables)

for config in confs:
    config.alltests.parent.mkdir(parents=True, exist_ok=True)
    with open(config.alltests, 'wb') as stream:
        pickle.dump(pooled[config.ind], stream)
