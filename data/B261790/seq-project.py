from collections import defaultdict
from pathlib import Path

from biobit.toolkit import nfcore, seqproj

FOLDER = Path(__file__).parent

FASTQ = FOLDER / "fastq"

loaded = nfcore.fetchngs.load_seqproj(
    samplesheet=FASTQ / "samplesheet.csv",
    fastq_root=Path("fastq")
)
assert len(loaded) == 1, f"Expected 1 seqproject, got {len(loaded)} for {FOLDER.name}"
project = loaded[0]

# Recreate the samples and experiments structure
samples = defaultdict(list)
for exp in project.experiments:
    mef, condition, selection, postfix = exp.sample.attributes.pop('title').split(', ')
    assert mef == "MEF"
    assert condition in {"mock", "HSV-1"}
    assert selection in {"IgG", "Z22", "IgG-input", "Z22-input", "RIP-input"}
    assert postfix.endswith(f" [{FOLDER.name}]")
    replica = postfix.split(' ')[0]
    assert replica in {"1", "2"}

    exp.sample.attributes["cells"] = "MEF"
    exp.sample.attributes["treatment"] = condition
    exp.sample.attributes["replica"] = replica

    ip = selection.replace("-input", "")
    samples[condition, ip, replica].append(exp)

    object.__setattr__(exp.library, "source", {"RNA"})
    object.__setattr__(exp.library, "strandedness", seqproj.Strandedness.Reverse)
    if selection in {"IgG", "Z22"}:
        object.__setattr__(exp.library, "selection", {"rRNA depletion", f"{selection} RIP"})
        exp.library.attributes["RIP"] = selection
    else:
        assert selection.endswith("-input")
        object.__setattr__(exp.library, "selection", {"rRNA depletion"})
        exp.library.attributes["RIP"] = "input"

    title = f"MEF_{condition}_{selection}_{replica}"
    exp.attributes['title'] = title

# Merge IP and input experiments into one biological sample to better reflect the experimental design
new_samples = []
for (condition, ip, replica), exps in samples.items():
    # Keep extra RIP inputs as is
    if ip == "RIP":
        assert len(exps) == 1
        new_samples.append(exps[0].sample)
        continue

    # Otherwise, we expect exactly one IP and one input experiment per sample
    assert len(exps) == 2, f"Expected 2 experiment per sample, got {len(exps)} for {(condition, ip, replica)}"
    ip, inp = exps[0].sample, exps[1].sample
    assert ip.attributes == inp.attributes
    assert ip.organism == inp.organism

    sample = seqproj.Sample(
        ind=f"{ip.ind}-{inp.ind}", organism=ip.organism, attributes=ip.attributes
    )
    object.__setattr__(exps[0], "sample", sample)
    object.__setattr__(exps[1], "sample", sample)
    new_samples.append(sample)

object.__setattr__(project, "samples", new_samples)
object.__setattr__(project, "ind", FOLDER.name)
object.__setattr__(project, "description", "Z22 RIP-seq of HSV-1 infected MEFs (batch 1)")

seqproj.adapter.yaml.dump(project, FOLDER / "seq-project.yaml")
