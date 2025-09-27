import pickle

import zhuntrs
from joblib import Parallel, delayed

import ld

PARALLEL = Parallel(n_jobs=-1, verbose=100)

# # The code below can be used to reproduce thresholds calculation
# import numpy as np
# import itertools
#
#
# def zhscore(seq, dnlim):
#     mindn = dnlim[0]
#     maxdn = min(len(seq) // 2, dnlim[1])
#     assert dnlim[0] <= mindn <= maxdn <= dnlim[1]
#     starts, ends, zhscores, _, _ = zhuntrs.predict(
#         ''.join(seq).encode("ASCII"), mindn=mindn, maxdn=maxdn, threshold=0, wrap=True
#     )
#     assert starts[0] == 0 and ends[0] <= len(seq), (starts, ends, zhscores)
#     return zhscores[0]
#
#
# # Build the scores distribution
# allscores = {}
#
# PARALLEL.batch_size = 1024
# for k in range(ld.MFE.dnlim[0] * 2, ld.MFE.dnlim[1] * 2 + 1):
#     print("Processing dn", k)
#     allscores[k] = PARALLEL(delayed(zhscore)(seq, ld.MFE.dnlim) for seq in itertools.product("ACGT", repeat=k))
#
# simplified = {}
# for k, scores in allscores.items():
#     total = len(scores)
#     print(f"[K={k}] Before removing identical scores: {len(scores)}")
#     scores = set(round(x, 3) for x in scores)
#     print(f"[K={k}] After removing identical scores: {len(scores)} ({len(scores) / total:.2%})")
#     simplified[k] = np.array(sorted(scores))
#
# # Calculate the thresholds
# quantiles = [0.9, 0.95, 0.99]
# allthresholds = {k: np.quantile(scores, quantiles).round(3) for k, scores in simplified.items()}
# for ind, q in enumerate(quantiles):
#     thresholds = {k: v[ind] for k, v in allthresholds.items()}
#     print(f"Quantile {q}:\n\tthresholds = {thresholds}")
#
# Quantile 0.9:
# 	thresholds = {6: 8.245, 7: 8.245, 8: 20.95, 9: 20.95, 10: 70.321, 11: 70.321, 12: 407.54}
# Quantile 0.95:
# 	thresholds = {6: 18.456, 7: 18.456, 8: 46.596, 9: 46.596, 10: 169.397, 11: 169.397, 12: 1161.065}
# Quantile 0.99:
# 	thresholds = {6: 122.484, 7: 122.484, 8: 403.258, 9: 403.258, 10: 1117.132, 11: 1117.132, 12: 7366.93}


THRESHOLDS = {
    6: 18.456, 7: 18.456, 8: 46.596, 9: 46.596, 10: 169.397, 11: 169.397, 12: 1161.065
}
assert ld.MFE.zh_quantile == 0.95


def job(fold, thresholds, saveto):
    with open(fold, 'rb') as stream:
        *meta, seq, fold, stems = pickle.load(stream)

    seq = seq.replace('U', 'T')
    assert set(seq) == {'A', 'C', 'G', 'T'}
    seq = seq.encode("ASCII")

    zstems, astems = [], []
    for left, right in stems:
        assert left.len() == right.len()
        if left.len() < ld.MFE.dnlim[0] * 2:
            astems.append((left, right))
            continue

        # Fetch the sequence and run the predictions
        zhseq = seq[left.start: left.end]
        dn = min(left.len() // 2, ld.MFE.dnlim[1])
        assert ld.MFE.dnlim[0] <= dn <= ld.MFE.dnlim[1]

        starts, ends, zhscores, windows, _ = zhuntrs.predict(
            zhseq, mindn=dn, maxdn=dn, threshold=0, wrap=True
        )

        # Select the max score for windows that include the stem only
        maxzh = 0
        for start, end, zhscore in zip(starts, ends, zhscores):
            if start >= 0 and end <= left.len():
                assert end - start == dn * 2
                maxzh = max(maxzh, zhscore)

        if maxzh >= thresholds[dn * 2]:
            zstems.append((left, right))
        else:
            astems.append((left, right))

    with open(saveto, 'wb') as stream:
        pickle.dump((*meta, seq, fold, astems, zstems), stream)


PARALLEL.batch_size = 1
PARALLEL(
    delayed(job)(fold, THRESHOLDS, ld.MFE.root / (fold.name.split('.')[0] + '.scored-fold.pkl'))
    for fold in ld.MFE.root.glob("*.raw-fold.pkl")
)

for file in ld.MFE.root.glob("*.scored-fold.pkl"):
    with open(file, 'rb') as stream:
        start, end, strand, seq, fold, stems, zstems = pickle.load(stream)
    print(file.name)
    print(f"{start}-{end}({strand})")
    print(seq)
    print(fold)
    print(f'Z-stems: {len(zstems)}')
    for left, right in zstems:
        print(f"\t {left.start}-{left.end}->{right.start}-{right.end}: "
              f"{seq[left.start: left.end]} -> {seq[right.start: right.end]}")
    print()
