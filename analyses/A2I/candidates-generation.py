import copy
import pickle
from collections import defaultdict
from itertools import chain
from pathlib import Path

from biobit.core import ngs
from biobit.toolkit import seqproj
from joblib import Parallel, delayed, cpu_count

import ld
from assemblies import GRCm39, CHM13v2


def parse_seqproj(experiment: seqproj.Experiment) -> tuple[str, str, bool]:
    bam = experiment.attributes["__nfcore_rnaseq_bam__"]
    layout = experiment.ngs()
    if isinstance(layout, ngs.Layout.Paired):
        pe = True
        if layout.strandedness == ngs.Strandedness.Forward:
            stranding = "s/f"
        elif layout.strandedness == ngs.Strandedness.Reverse:
            stranding = "f/s"
        else:
            raise ValueError(f"Unstranded library: {layout.strandedness}")
    elif isinstance(layout, ngs.Layout.Single):
        pe = False
        if layout.strandedness == ngs.Strandedness.Forward:
            stranding = "s"
        elif layout.strandedness == ngs.Strandedness.Reverse:
            stranding = "f"
        else:
            raise ValueError(f"Unstranded library: {layout.strandedness}")
    else:
        raise ValueError(f"Unknown layout: {layout}")
    return bam, stranding, pe


workload = []
all_experiments = defaultdict(list)
for assembly in CHM13v2, GRCm39:
    for project in ld.c.reat.experiments[assembly.name]:
        prjind = project.ind.replace(" ", "_").replace("/", "_")
        pooled = defaultdict(list)
        for experiment in project.experiments:
            expind = experiment.ind.replace(" ", "_").replace("/", "_")

            # Run for each experiment
            saveto = ld.c.reat.candidates / assembly.name / prjind / f"{expind}.csv"
            saveto.parent.mkdir(parents=True, exist_ok=True)

            bam, stranding, paired = parse_seqproj(experiment)
            workload.append((assembly.fasta, [bam], stranding, paired, saveto, ld.c.VIRAL))

            # Derive the pooled key
            key = sorted(chain(
                experiment.attributes.items(), experiment.sample.attributes.items(),
                experiment.library.attributes.items(),
            )) + sorted(experiment.library.selection)
            key = tuple(x for x in key if x[0] not in {"replica", "title"} and not x[0].startswith("__"))

            # Record the experiment
            all_experiments[project.ind].append(experiment)
            pooled[key].append(experiment)

        # Add pooled samples
        for key, experiments in pooled.items():
            assert len(experiments) > 1, key

            # Sample is a deepcopy with new index and "pooled" replica
            sample = copy.deepcopy(experiments[0].sample)
            object.__setattr__(sample, "ind", f"+".join(x.sample.ind for x in experiments) + "[merged]")
            sample.attributes["replica"] = "pooled"

            # Library is just a deepcopy
            library = copy.deepcopy(experiments[0].library)

            # Experiment has no fields except required
            runs = tuple(chain(*[x.runs for x in experiments]))
            experiment = seqproj.Experiment("+".join(x.ind for x in experiments) + "[merged]", sample, library, runs)

            # Add the pooled sample to the workload
            bam, stranding, paired = parse_seqproj(experiments[0])
            bams = [bam]
            for exp in experiments[1:]:
                bam, _stranding, _paired = parse_seqproj(exp)
                bams.append(bam)
                assert stranding == _stranding and paired == _paired

            expind = experiment.ind.replace(" ", "_").replace("/", "_")
            saveto = ld.c.reat.candidates / assembly.name / prjind / f"{expind}.csv"
            workload.append((assembly.fasta, bams, stranding, paired, saveto, ld.c.VIRAL))

            # Record the pooled experiment
            all_experiments[project.ind].append(experiment)

ld.c.reat.seqproj.parent.mkdir(parents=True, exist_ok=True)
with open(ld.c.reat.seqproj, 'wb') as stream:
    pickle.dump(dict(all_experiments), stream)

filtered = []
for w in workload:
    if all(Path(x).exists() for x in w[1]):
        filtered.append(w)
    else:
        print(f"Skipped: ", w)
workload = filtered

# Run the initial REAT calculation
Parallel(n_jobs=2)(delayed(ld.utils.reat.run)(*args, threads=cpu_count()) for args in workload)
