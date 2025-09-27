import pickle
from collections import defaultdict

from stories.RIP import pcalling
from stories.RIP.clustering import ld


def resolve(name: str, comparisons: list[pcalling.Config]) -> ld.Config:
    assert name in ("Z-RNA[CHM13v2]", "Z-RNA[GRCm39]"), name
    strong_min_replication = {"Z-RNA[CHM13v2]": 8, "Z-RNA[GRCm39]": 4}[name]
    peaks = ld.PeaksConfig(
        min_length=36, coverage_gaps_tolerance=5_000,
        strong_min_replication=strong_min_replication, relaxed_min_replication=2,
        relaxed_min_editing_sites=2
    )

    min_replication = {"Z-RNA[CHM13v2]": 2, "Z-RNA[GRCm39]": 2}[name]
    dsRNA = ld.dsRNAConfig(min_replication=min_replication)
    groups = ld.ClusteringConfig()

    return ld.Config(
        name, tuple(comparisons), ld.RESULTS / name, peaks, dsRNA, groups
    )


pcallconf = pcalling.Config.load()

# Group configs by the target
groups = defaultdict(list)
for conf in pcallconf:
    if len(conf.signal) > 1 or len(conf.control) > 1:
        print(f"Skipping {conf.ind} due to multiple signal "
              f"(N={len(conf.signal)})/control (N={len(conf.control)}) libraries")
        continue

    if "MEF" in conf.ind:
        groups["Z-RNA[GRCm39]"].append(conf)
    else:
        assert "HT-29" in conf.ind
        groups["Z-RNA[CHM13v2]"].append(conf)

allconfigs = []
for target, configs in groups.items():
    allconfigs.append(resolve(target, configs))

with open(ld.config.PKL, 'wb') as stream:
    pickle.dump(allconfigs, stream)
