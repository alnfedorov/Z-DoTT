import pickle
from typing import Iterable

import pandas as pd
from biobit.toolkit import seqproj

import ld
from stories.aberrantome import Config


def resolve(assembly: str, prj: seqproj.Project) -> Iterable[Config]:
    root = ld.RESULTS / assembly
    df = pd.DataFrame(
        [{"Experiment": exp, **exp.sample.attributes, **exp.library.attributes} for exp in prj.experiments]
    )
    match prj.ind:
        case 'HT-29 HSV-1/IAV [B831009, JCC411]':
            results = root / "B831009"

            # Drop IPs
            mask = df['Experiment'].apply(lambda x: x.library.selection == {"Total RNA", "rRNA depletion"})
            df = df[mask]
            control = df.loc[df['treatment'] == 'mock', 'Experiment'].to_list()
            assert len(control) == 4

            for virus in "HSV-1", "IAV":
                treatment = df.loc[df['treatment'] == virus, 'Experiment'].to_list()
                assert len(treatment) == 8
                yield Config(
                    f"HT-29: {virus} vs mock", assembly, tuple([(prj.ind, exp) for exp in treatment]),
                    tuple([(prj.ind, exp) for exp in control]), virus, "Mock", results
                )
        case 'PRJNA256013':
            results = root / "PRJNA256013"
            for protocol in "Total-RNA", "4sU-RNA":
                subdf = df[df['protocol'] == protocol]

                control = subdf.loc[subdf['time-point'] == 'mock', 'Experiment'].to_list()
                assert len(control) == 2

                for time_point in subdf['time-point'].unique():
                    if time_point == 'mock':
                        continue

                    treatment = subdf.loc[subdf['time-point'] == time_point, 'Experiment'].to_list()
                    assert len(treatment) == 2
                    yield Config(
                        f"HFF: {protocol} {time_point} vs mock", assembly, tuple([(prj.ind, exp) for exp in treatment]),
                        tuple([(prj.ind, exp) for exp in control]), time_point, "Mock", results
                    )
        case 'MEF IAV batch 1 [B256178, JCC400]':
            results = root / "B256178"

            # Drop IPs
            mask = df['Experiment'].apply(lambda x: x.library.selection == {"Total RNA", "rRNA depletion"})
            df = df[mask]
            control = df.loc[df['treatment'] == 'mock', 'Experiment'].to_list()
            assert len(control) == 6

            treatment = df.loc[df['treatment'] == 'IAV', 'Experiment'].to_list()
            assert len(treatment) == 6
            yield Config(
                "MEF: IAV vs mock", assembly, tuple([(prj.ind, exp) for exp in treatment]),
                tuple([(prj.ind, exp) for exp in control]), "IAV", "Mock", results
            )
        case 'MEF HSV-1 batch 1 [B261790, JCC411]':
            result = root / "B261790"

            # Drop IPs
            mask = df['Experiment'].apply(lambda x: x.library.selection == {"Total RNA", "rRNA depletion"})
            df = df[mask]
            control = df.loc[df['treatment'] == 'mock', 'Experiment'].to_list()
            assert len(control) == 6

            treatment = df.loc[df['treatment'] == 'HSV-1', 'Experiment'].to_list()
            assert len(treatment) == 6

            yield Config(
                "MEF: HSV-1 vs mock", assembly, tuple([(prj.ind, exp) for exp in treatment]),
                tuple([(prj.ind, exp) for exp in control]), "HSV-1", "Mock", result
            )
        case 'PRJEB75711':
            result = root / "PRJEB75711"

            for (fraction, time_point), subdf in df.groupby(['fraction', 'time-point']):
                control = subdf.loc[subdf['IAV'] == 'mock', 'Experiment'].to_list()
                assert len(control) == 3, (fraction, time_point, control)

                for IAV in "WT", "dNS1":
                    treatment = subdf.loc[subdf['IAV'] == f'IAV-{IAV}', 'Experiment'].to_list()
                    assert len(treatment) == 3, (fraction, time_point, IAV, treatment)

                    yield Config(
                        f"A549: IAV {IAV}-{time_point} {fraction}", assembly,
                        [(prj.ind, exp) for exp in treatment], [(prj.ind, exp) for exp in control],
                        f"IAV-{IAV} {time_point} [{fraction}]", f"Mock [{fraction}]", result
                    )
        case 'PRJNA637636':
            result = root / "PRJNA637636"

            control = df.loc[df['HSV-1'] == 'mock', 'Experiment'].to_list()
            assert len(control) == 2

            for HSV1 in df['HSV-1'].unique():
                if HSV1 == 'mock' or HSV1 == 'delta ICP4':
                    continue

                treatment = df.loc[df['HSV-1'] == HSV1, 'Experiment'].to_list()
                assert len(treatment) == 2, (HSV1, treatment)

                yield Config(
                    f"HFF: HSV-1 {HSV1} vs mock", assembly, [(prj.ind, exp) for exp in treatment],
                    [(prj.ind, exp) for exp in control], f"HSV-1 {HSV1}", "Mock", result
                )
        case _:
            raise ValueError(prj.ind)


configs = []
for assembly, projects in ld.comparisons.series.items():
    for prj in projects:
        configs.extend(resolve(assembly, prj))

# Short report
for cmp in configs:
    print(cmp.ind)
    print(f"\t{cmp.assembly}")
    print("\tSignal:")
    for _, trt in cmp.treatment:
        print("\t\t", trt.attributes['title'])
    print("\tControl:")
    for _, ctr in cmp.control:
        print("\t\t", ctr.attributes['title'])
    print()

ld.comparisons.pkl.parent.mkdir(parents=True, exist_ok=True)
with open(ld.comparisons.pkl, 'wb') as f:
    pickle.dump(configs, f)
