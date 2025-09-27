import pickle
from collections import defaultdict
from typing import Any

from biobit.toolkit import seqproj, reaper

from stories.RIP.pcalling.ld import config, RESULTS
from stories.nextflow.series import internal as nextflow


def resolve(
        prj: seqproj.Project
) -> tuple[list[config.Config], dict[Any, list[config.Config]]]:
    prjind = {
        "HT-29 HSV-1/IAV [B831009, JCC411]": "B831009",
        "MEF IAV batch 1 [B256178, JCC400]": "B256178",
        "MEF HSV-1 batch 1 [B261790, JCC411]": "B261790",
        "MEF HSV-1 batch 2 [B319096, JCC478]": "B319096",
    }[prj.ind]

    # Group all things by sample
    samples: dict[str, list[seqproj.Experiment]] = defaultdict(list)
    for exp in prj.experiments:
        samples[exp.sample.ind].append(exp)

    root = RESULTS / prjind

    model = reaper.model.RNAPileup() \
        .set_min_signal(10) \
        .set_sensitivity(1e-6)

    cutoff = reaper.pcalling.ByCutoff() \
        .set_cutoff(2.0) \
        .set_min_length(6) \
        .set_merge_within(6)

    nms = reaper.postfilter.NMS() \
        .set_group_within(10_000) \
        .set_fecutoff(2.0) \
        .set_slopfrac(1.0) \
        .set_sloplim(1_000, 10_000) \
        .set_sensitivity(1e-6)

    # Derive individual and/or pooled comparisons
    comparisons, pooled = [], defaultdict(list)
    for experiments in samples.values():
        # Derive the selection for each experiment from the given sample
        by_selection = {}
        for exp in experiments:
            selection = exp.library.selection - {'Total RNA', 'rRNA depletion'}
            assert len(selection) <= 1, selection
            selection = selection.pop() if selection else 'Input'
            assert selection not in by_selection, (prj.ind, by_selection, selection)
            by_selection[selection] = exp

        # Comparison: matched IP & input
        control = by_selection.pop("Input")
        organism = frozenset(control.sample.organism)
        host = "Homo sapiens" if "Homo sapiens" in organism else "Mus musculus"
        assembly = {"Homo sapiens": "CHM13v2", "Mus musculus": "GRCm39"}[host]
        for Ab in "FLAG RIP", "Z22 RIP",:  # "IgG RIP":
            if Ab not in by_selection:
                continue
            signal = by_selection.pop(Ab)
            title = f"{prjind}+{signal.attributes['title']} vs {control.attributes['title']}"

            if Ab == "FLAG RIP" and "mock" in title:
                continue

            comparisons.append(config.Config(
                title, prj.ind, host, organism, assembly, (signal,), (control,), model, cutoff, nms, root
            ))
            pooled[signal.attributes['title'][:-2]].append(comparisons[-1])

    # for title, cmps in pooled.items():
    #     assert 2 <= len(cmps) <= 4, title
    #     signal = tuple(chain(*(cmp.signal for cmp in cmps)))
    #     control = tuple(chain(*(cmp.control for cmp in cmps)))
    #
    #     host, organism, assembly = cmps[0].host, cmps[0].organism, cmps[0].assembly
    #     assert all(x.host == host and x.organism == organism and x.assembly == assembly for x in cmps)
    #     comparisons.append(config.Config(
    #         title, prj.ind, host, organism, assembly, signal, control, model, cutoff, nms, root
    #     ))

    # Derive the group identifier
    match prjind:
        case "B831009":
            group = "Z-RNA[CHM13v2]"
        case "B319096" | "B261790" | "B256178":
            group = "Z-RNA[GRCm39]"
        case _:
            raise ValueError(f"Unknown project ind: {prjind}")

    return comparisons, {group: comparisons}


# Resolve all pairwise comparisons
comparisons = []
for prj in [
    nextflow.B256178,  # IAV infection of MEFs (batch 1, Z22 + FLAG RIP-seq)
    nextflow.B261790,  # HSV-1 infection of MEFs (batch 1, Z22 RIP-seq)
    nextflow.B319096,  # HSV-1 infection of MEFs (batch 2, FLAG RIP-seq)
    nextflow.B831009  # HSV-1/IAV RIP from HT-29s (Z22 + FLAG RIP-seq)
]:
    cmp, grouped = resolve(prj)
    comparisons.extend(cmp)

# Short report
for cmp in comparisons:
    print(cmp.ind)
    print(f"\t {cmp.project}: {cmp.host} ({tuple(cmp.organism)}")
    print("\tSignal:")
    for signal in cmp.signal:
        print("\t\t", signal.attributes['title'])
    print("\tControl:")
    for control in cmp.control:
        print("\t\t", control.attributes['title'])
    print()

with open(config.PKL, 'wb') as f:
    pickle.dump(comparisons, f, protocol=pickle.HIGHEST_PROTOCOL)
