import pickle

import pandas as pd
from joblib import Parallel, delayed

import ld
from assemblies import CHM13v2, GRCm39
from stories.aberrantome import Config

ANNOTATION = {
    "CHM13v2": CHM13v2.gencode.load(),
    "GRCm39": GRCm39.gencode.load()
}


def job(cmp: Config):
    annotation = ANNOTATION[cmp.assembly]

    for fname in cmp.read_through, cmp.read_in, cmp.divergent:
        df = pd.read_csv(fname)
        df['Name'] = df['Ensembl ID'].apply(lambda ind: annotation.genes[ind].attrs.name)
        df.to_csv(fname, index=False)

    for fname in cmp.intronic,:
        df = pd.read_csv(fname)

        names = []
        for eid in df['Ensembl ID']:
            eid, rank = eid.split('-', maxsplit=1)
            names.append(
                f"{annotation.genes[eid].attrs.name}[#{rank}]"
            )
        df['Name'] = names
        df.to_csv(fname, index=False)


with open(ld.comparisons.pkl, 'rb') as stream:
    comparisons = pickle.load(stream)

Parallel(n_jobs=-1, backend='multiprocessing')(delayed(job)(cmp) for cmp in comparisons)
