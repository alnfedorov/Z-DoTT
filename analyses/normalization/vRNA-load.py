import pandas as pd

import ld
from stories import normalization

pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

counts = normalization.load_counts()

RESULTS = []

for assembly, cnts in counts.items():
    projects = normalization.SERIES[assembly]
    items = {(prj.ind, exp.ind): exp for prj in projects for exp in prj.experiments}

    virus = ['HSV-1[+]', 'HSV-1[-]', 'IAV[+]', 'IAV[-]']
    host = [x for x in cnts.columns if any(x.startswith(prefix) for prefix in ('ENS', 'chr', 'GL', 'JH', 'MU'))]
    assert set(cnts.columns) == set(host + virus + ['source']), \
        f"Unexpected columns: {set(cnts.columns) - set(host + virus + ['source'])}"

    cnts['Host [reads]'] = cnts[host].sum(axis=1)
    cnts['Virus [reads]'] = cnts[virus].sum(axis=1)

    total = cnts['Host [reads]'] + cnts['Virus [reads]']
    cnts['Host [proportion]'] = cnts['Host [reads]'] / total
    cnts['Virus [proportion]'] = cnts['Virus [reads]'] / total

    cnts = cnts.drop(columns=virus + host)
    RESULTS.append(cnts)

RESULTS = pd.concat(RESULTS, axis=0)
RESULTS.to_pickle(ld.VIRAL_RNA_LOAD)
