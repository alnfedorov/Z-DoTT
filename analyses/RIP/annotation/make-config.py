import pickle
import warnings
from collections.abc import Iterator

import pandas as pd

import ld
from stories.RIP import annotation, clustering


def resolve(clconf: clustering.Config) -> Iterator[annotation.Config]:
    # Load all partitions
    if clconf.clusters.results.exists():
        with open(clconf.clusters.results, "rb") as stream:
            elements = pickle.load(stream)
    else:
        warnings.warn(f"No groups found for {clconf.ind}")
        elements = []

    # Derive tests for each original group
    comparisons = [cmp for cmp in clconf.comparisons if len(cmp.signal) == 1]  # Ignore pooled comparisons
    match clconf.ind:
        case "Z-RNA[CHM13v2]":
            features = []
            for cmp in comparisons:
                for exp in cmp.signal + cmp.control:
                    record = {'sample': exp.sample.ind, "treatment": exp.sample.attributes["treatment"]}
                    RIP = exp.library.selection - {'Total RNA', 'rRNA depletion'}
                    assert len(RIP) <= 1, RIP
                    record['RIP'] = RIP.pop() if RIP else 'input'
                    record['ind'] = (cmp.project, exp.ind)
                    features.append(record)

            features = pd.DataFrame(features).set_index('ind').drop_duplicates()
            features["RIP"] = features["RIP"].map({"Z22 RIP": "Z22", "input": "input", "FLAG RIP": "FLAG"})
            features["treatment"] = features["treatment"].map({"mock": "mock", "HSV-1": "HSV1", "IAV": "IAV"})
            features["group"] = [f"{rip}.{treatment}" for rip, treatment in zip(features['RIP'], features['treatment'])]

            tests = []
            for baseline, names in [
                ("input.mock", {
                    "Z22 enriched [mock]": "group_Z22.mock_vs_input.mock",
                    "Input induced [HSV-1]": "group_input.HSV1_vs_input.mock",
                    "Input induced [IAV]": "group_input.IAV_vs_input.mock",
                }),
                ("Z22.mock", {
                    "Z22 induced [HSV-1]": "group_Z22.HSV1_vs_Z22.mock",
                    "Z22 induced [IAV]": "group_Z22.IAV_vs_Z22.mock",
                }),
                ("input.HSV1", {
                    "Z22 enriched [HSV-1]": "group_Z22.HSV1_vs_input.HSV1",
                    "FLAG enriched [HSV-1]": "group_FLAG.HSV1_vs_input.HSV1",
                }),
                ("input.IAV", {
                    "Z22 enriched [IAV]": "group_Z22.IAV_vs_input.IAV",
                    "FLAG enriched [IAV]": "group_FLAG.IAV_vs_input.IAV",
                })
            ]:
                for name, cmp in names.items():
                    alternative, alpha = {
                        "induced": ("greaterAbs", 0.1),
                        "enriched": ("greater", 0.01)
                    }[name.split()[1]]
                    tests.append(annotation.StatTest(
                        features, "~group", alternative, {"group": baseline}, {name: cmp}, log2fc=0.0, alpha=alpha
                    ))
        case "Z-RNA[GRCm39]":
            features = []
            for cmp in comparisons:
                for exp in cmp.signal + cmp.control:
                    record = {'sample': exp.sample.ind, "treatment": exp.sample.attributes["treatment"]}
                    RIP = exp.library.selection - {'Total RNA', 'rRNA depletion'}
                    assert len(RIP) <= 1, RIP
                    record['RIP'] = RIP.pop() if RIP else 'input'
                    record['project'] = cmp.project
                    record['ind'] = (cmp.project, exp.ind)
                    features.append(record)

            features = pd.DataFrame(features).set_index('ind').drop_duplicates()
            features["RIP"] = features["RIP"].apply(lambda x: x.replace(" RIP", ""))
            features["treatment"] = features["treatment"].map({"mock": "mock", "HSV-1": "HSV1", "IAV": "IAV"})
            features['group'] = [f"{rip}.{treatment}" for rip, treatment in zip(features['RIP'], features['treatment'])]
            features['batch'] = ["_".join(x.split()[:4]) for x in features['project']]

            # Renumerate samples inside each project
            tests = []
            for baseline, names in [
                ('input.mock', {
                    'Z22 enriched [mock]': 'group_Z22.mock_vs_input.mock',
                    'Input induced [IAV]': 'group_input.IAV_vs_input.mock',
                    'Input induced [HSV-1]': 'group_input.HSV1_vs_input.mock',
                }),
                ('Z22.mock', {
                    'Z22 induced [IAV]': 'group_Z22.IAV_vs_Z22.mock',
                    'Z22 induced [HSV-1]': 'group_Z22.HSV1_vs_Z22.mock',
                }),
                ('input.IAV', {
                    'Z22 enriched [IAV]': 'group_Z22.IAV_vs_input.IAV',
                    'FLAG enriched [IAV]': 'group_FLAG.IAV_vs_input.IAV',
                }),
                ('input.HSV1', {
                    'Z22 enriched [HSV-1]': 'group_Z22.HSV1_vs_input.HSV1',
                    'FLAG enriched [HSV-1]': 'group_FLAG.HSV1_vs_input.HSV1',
                }),
            ]:
                for name, cmp in names.items():
                    alternative, alpha = {
                        "induced": ("greaterAbs", 0.1),
                        "enriched": ("greater", 0.01)
                    }[name.split()[1]]
                    tests.append(annotation.StatTest(
                        features, "~batch + group", alternative, {"group": baseline}, {name: cmp},
                        log2fc=0.0, alpha=alpha
                    ))
        case _:
            raise ValueError(f"Unknown grouping: {clconf.ind}")

    yield annotation.Config(
        clconf.ind, clconf.ind, tuple(elements), tuple(comparisons), tuple(tests), ld.RESULTS / clconf.ind
    )


configs = []
for conf in clustering.Config.load():
    configs.extend(resolve(conf))

with open(ld.config.PKL, "wb") as stream:
    pickle.dump(configs, stream)
