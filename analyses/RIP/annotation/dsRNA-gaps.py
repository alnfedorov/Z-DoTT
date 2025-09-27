import pickle
from collections import defaultdict

from joblib import Parallel, delayed
from tqdm import tqdm

import ld
from stories.RIP import annotation, pcalling

SIGNAL = ld.cache.signal.load()


def job(config: annotation.Config, cmp: pcalling.Config):
    signal = SIGNAL[config.ind]

    results = {}
    for partition in tqdm(config.elements):
        cache = signal[partition.contig, partition.orientation][cmp.project, cmp.ind]

        # Weight recorded inner_gaps by the IP signal
        scores = {}
        for rep in partition.invrep:
            distance = rep.inner_gap()

            maxsignal = float('-inf')
            for arm in rep.segments:
                maxsignal = max(maxsignal, min(cache[arm.left], cache[arm.right]))

            if distance in scores:
                scores[distance] = max(maxsignal, scores[distance])
            else:
                scores[distance] = maxsignal

        # Use the distance with the highest signal
        if scores:
            results[partition.ind] = max(scores.items(), key=lambda x: x[1])[0]
        else:
            results[partition.ind] = None
    return config.ind, cmp.ind, results


configs = annotation.Config.load()
_results = Parallel(n_jobs=-1, verbose=100, backend='sequential')(
    delayed(job)(config, cmp) for config in configs for cmp in config.comparisons
)

# Group the results and save to the cache
results = defaultdict(dict)
for conf, cmp, res in _results:
    results[conf][cmp] = res

for conf in configs:
    conf.loops.parent.mkdir(parents=True, exist_ok=True)
    with open(conf.loops, 'wb') as stream:
        pickle.dump(results[conf.ind], stream)
