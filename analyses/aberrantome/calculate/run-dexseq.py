from itertools import chain
from pathlib import Path
from subprocess import check_call
from tempfile import NamedTemporaryFile

import pandas as pd
from joblib import Parallel, delayed
from statsmodels.stats.multitest import fdrcorrection

import ld
from stories import normalization
from stories.aberrantome import Config

DEXSEQ = Path(__file__).with_suffix('.R')

# Libraries normalization strategy
NORMALIZATION = normalization.median_of_ratios()

# Raw counts for each aberrant transcription region
COUNTS = {
    "read-through": pd.read_pickle(ld.counts.read_through),
    "read-in": pd.read_pickle(ld.counts.read_in),
    "divergent": pd.read_pickle(ld.counts.divergent),
    "intronic": pd.read_pickle(ld.counts.intronic)
}

# Calculate FPKMs based on host-mapped reads only
TOTAL_READS = normalization.load_vRNA_load().set_index('source')['Host [reads]'].to_dict()


def job(cmp: Config, cmpkey: str):
    trtkey = list(cmp.treatment_key())
    ctrlkey = list(cmp.control_key())
    allkeys = trtkey + ctrlkey

    # Calculate FPKMs for each region
    fpkm = COUNTS[cmpkey][allkeys + ['Ensembl ID', 'Length', 'Region']].copy()
    for column in allkeys:
        fpkm[column] = (fpkm[column] * 1e9) / (fpkm['Length'] * TOTAL_READS[column])

    # Create a list of reference regions passing the threshold in each sample
    references = fpkm.loc[fpkm['Region'] == 'reference', ['Ensembl ID'] + allkeys].set_index('Ensembl ID')
    passed = (references > ld.thr.scores.min_ref_fpkm).sum(axis=1)
    passed = set(references[passed >= len(allkeys)].index)

    print(f"Passed: {len(passed)} ({len(passed) / len(references):.2%}) of all reference regions")

    # Filter out regions with insufficient reference coverage
    fpkm = fpkm[fpkm['Ensembl ID'].isin(passed)].copy()

    # Calculate aberrant transcription scores
    scores = pd.pivot_table(fpkm, index='Ensembl ID', columns='Region', values=allkeys)
    for key in allkeys:
        aberrant, reference = scores[(key, 'aberrant')], scores[(key, 'reference')]
        scores[(key, 'score')] = aberrant / (aberrant + reference)

    for group, name in zip([trtkey, ctrlkey], ['Score [trt]', 'Score [ctrl]']):
        group = [(x, 'score') for x in group]
        scores[name] = scores[group].median(axis=1, skipna=True)

    columns = ['Ensembl ID', 'Score [trt]', 'Score [ctrl]']
    scores = scores.reset_index()[columns].dropna(how='any')
    scores.columns = columns

    # Fetch counts
    counts = COUNTS[cmpkey][allkeys + ['Ensembl ID', 'Region']]
    counts = counts[counts['Ensembl ID'].isin(passed)].copy()

    # Normalize counts
    _, _, scaling = NORMALIZATION[cmp.assembly].scaling_factors(
        {(prj, exp.ind): [(prj, exp.ind)] for prj, exp in chain(cmp.treatment, cmp.control)},
    )
    for key, val in scaling.to_dict().items():
        counts[key] = (counts[key] / val).round(0).astype(int)

    # Create the samples table
    allsamples, rename = [], {}
    for samples, condition in (trtkey, 'treatment'), (ctrlkey, 'control'):
        for ind, key in enumerate(samples):
            rename[key] = f"{condition}_{ind}"
            allsamples.append({"sample": rename[key], "condition": condition})
    allsamples = pd.DataFrame(allsamples)
    counts = counts.rename(columns=rename)

    # Setup and run the DEXSeq
    with (
        NamedTemporaryFile() as smpltable,
        NamedTemporaryFile() as cntstable,
        NamedTemporaryFile() as saveto
    ):
        allsamples.to_csv(smpltable, index=False)

        counts = counts.rename(columns={"Ensembl ID": "group", "Region": "element"})
        counts.to_csv(cntstable, index=False)

        parameters = [
            "Rscript", DEXSEQ, smpltable.name, cntstable.name, saveto.name, 'control', str(True)
        ]
        print(parameters, "\n")
        check_call(parameters)
        results = pd.read_csv(saveto, compression="gzip", index_col=0)

    results = results.rename(columns={
        f"log2fold_treatment_control": "log2fold_change",
        "groupID": "Ensembl ID",
        "featureID": "Region",
        "exonBaseMean": "base_mean"
    }).drop(columns=['treatment', 'control']).reset_index(drop=True)

    # Select only aberrant regions, drop untested records, and recalculate padj
    results = results[results['Region'] == 'aberrant'].drop(columns='Region').dropna(how='any')
    results['padj'] = fdrcorrection(results['pvalue'])[1]

    # Merge scores and save results
    saveto = {
        "read-through": cmp.read_through,
        "read-in": cmp.read_in,
        "intronic": cmp.intronic,
        "divergent": cmp.divergent
    }[cmpkey]

    results = pd.merge(results, scores, on='Ensembl ID', how='inner')
    saveto.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(saveto, index=False)


comparisons = ld.comparisons.load()
Parallel(n_jobs=-1, backend='multiprocessing')(delayed(job)(cmp, key) for cmp in comparisons for key in COUNTS.keys())
