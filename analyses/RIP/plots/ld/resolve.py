from collections import defaultdict

import numpy as np
import pandas as pd
from biobit.deprecated.repmasker import RepmaskerClassification


def sequence(
        repcls: RepmaskerClassification,
        weights: dict[tuple[str, str] | tuple[tuple[str, str], tuple[str, str]], float]
) -> tuple[str, str]:
    assert all(v > 0 for v in weights.values()), weights
    assert abs(sum(weights.values()) - 1) < 1e-3, weights

    resolved = defaultdict(float)
    for k, v in weights.items():
        if isinstance(k[0], str):
            k = (k,)

        if len(k) == 1:
            key = ('Unresolved', 'Unresolved')
        else:
            assert len(k) == 2, k
            left, right = k
            if left[0] == right[0] == 'Repeat-free':
                key = ('Inv. non-repetitive', 'Inv. non-repetitive')
            elif left[0] == 'Repeat-free' or right[0] == 'Repeat-free':
                key = ('Complex', 'Complex')
            else:
                lname, lstrand = left
                rname, rstrand = right

                _, lfamily, lcls = repcls.classify(lname)
                _, rfamily, rcls = repcls.classify(rname)
                lfamily, rfamily = (lfamily if lfamily else lcls), (rfamily if rfamily else rcls)

                if lcls == rcls and lstrand != rstrand:
                    cls = f"Inv. {lcls}"
                elif lcls == rcls:
                    cls = lcls
                else:
                    cls = "Complex"

                if lfamily == rfamily and lstrand != rstrand:
                    assert cls == f"Inv. {lcls}"
                    family = f"Inv. {lfamily}"
                elif lfamily == rfamily:
                    assert cls == lcls
                    family = lfamily
                else:
                    family = "Complex"
                key = (cls, family)

        resolved[key] = max(resolved[key], v)

    resolved = sorted(resolved.items(), key=lambda x: (x[1], x[0]))
    solution, weight = resolved.pop()

    # If the solution is not related to inverted sequences, check if there is a similarly weighted inverted sequence
    if not solution[0].startswith("Inv.") or not solution[1].startswith("Inv."):
        while resolved:
            nxt, nxtweight = resolved.pop()
            if nxt[0].startswith("Inv.") and nxt[1].startswith('Inv.') and abs(nxtweight - weight) < 1e-2:
                solution = nxt
                break
    return solution


def by_refexp(df: pd.DataFrame, refexp: list[str], repcls: RepmaskerClassification) -> pd.DataFrame:
    # Take the most frequent location
    for postfix in "auto", "manual":
        locs = [f"location-{postfix} [{x}]" for x in refexp]
        loc = df[locs].apply(lambda x: x.value_counts().idxmax(), axis=1).to_list()
        df[f'Ensembl ID [{postfix}]'], df[f'Gene name [{postfix}]'], df[f'location [{postfix}]'] = zip(*loc)

    # Take the most frequent loop-size
    loop = [f"loop-size [{x}]" for x in refexp]
    df['loop-size'] = df[loop].fillna(-1).apply(lambda x: x.value_counts().idxmax(), axis=1)

    # Resolve the sequence for each experiment
    seqs = [f"sequence [{x}]" for x in refexp]
    for col in seqs:
        df[f'sequence resolved {col}'] = df[col].apply(lambda x: sequence(repcls, x))
    cols = [f'sequence resolved {col}' for col in seqs]
    df['sequence resolved'] = df[cols].apply(lambda x: x.value_counts().idxmax(), axis=1)

    df['sequence[cls]'] = df['sequence resolved'].apply(lambda x: x[0])
    df['sequence[family]'] = df['sequence resolved'].apply(lambda x: x[1])

    df = df.drop(columns=cols + ['sequence resolved'])
    return df


def virus_induction(df: pd.DataFrame, experiments: dict[str, list[str]], log2fc: float, padj: float) -> pd.DataFrame:
    for key, exps in experiments.items():
        mask = np.zeros(len(df), dtype=bool)
        for exp in exps:
            mask |= (df[(exp, 'log2FoldChange')] >= log2fc) & (df[(exp, 'padj')] <= padj)
        df[f'{key} induced'] = mask
    return df


def rip_enrichment(df: pd.DataFrame, experiments: dict[str, str], log2fc: float, padj: float) -> pd.DataFrame:
    rename = {}
    for newname, oldname in experiments.items():
        for col in "log2FoldChange", "padj", "pvalue":
            assert (oldname, col) not in rename, (newname, oldname, col)
            rename[oldname, col] = (newname, col)

        mask = (df[(oldname, 'log2FoldChange')] >= log2fc) & (df[(oldname, 'padj')] <= padj)
        df[newname] = mask
    df = df.rename(columns=rename)
    return df
