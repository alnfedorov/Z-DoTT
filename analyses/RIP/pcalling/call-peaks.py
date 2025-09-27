import copy
import pickle
from multiprocessing import Pool
from multiprocessing import cpu_count
from typing import Any

from biobit import io
from biobit.core.loc import Strand
from biobit.core.ngs import Layout
from biobit.toolkit import reaper as rp, nfcore

import ld
import utils
from stories import normalization

SAVERS = Pool(cpu_count())


def reader(path: str, layout: Layout) -> io.bam.Reader:
    assert isinstance(layout, Layout.Paired), f"Unsupported layout: {layout} ({path})"
    return io.bam.Reader(path, inflags=3, exflags=2572, minmapq=0)


def save_results(cmp: ld.Config, harvest: list[rp.HarvestRegion]) -> list[tuple[str, Any]]:
    total_peaks = sum(len(x.filtered_peaks) for x in harvest)
    print(f"[{cmp.ind}] Finished {cmp.ind} -> {total_peaks:,} peaks")

    # Construct & save the BED files
    paths = cmp.reaper
    for p in paths:
        p.parent.mkdir(parents=True, exist_ok=True)

    handles = []

    # Save segments
    for get_segments, saveto in [
        (lambda x: x.control, paths.control), (lambda x: x.signal, paths.signal), (lambda x: x.modeled, paths.modeled)
    ]:
        lines = []
        for region in harvest:
            contig = region.contig
            strand = str(region.orientation)
            for segment in get_segments(region):
                lines.append(f"{contig}\t{segment.start}\t{segment.end}\t.\t.\t{strand}\n")
        handles.append((saveto, SAVERS.apply_async(ld.write_gz, (saveto, lines))))

    # Save peaks
    for get_peaks, saveto in [
        (lambda x: x.raw_peaks, paths.raw_peaks), (lambda x: x.filtered_peaks, paths.filtered_peaks)
    ]:
        lines = []
        for region in harvest:
            contig = region.contig
            strand = str(region.orientation)
            for peak in get_peaks(region):
                lines.append(f"{contig}\t{peak.interval.start}\t{peak.interval.end}\t.\t{peak.value}\t{strand}\n")
        handles.append((saveto, SAVERS.apply_async(ld.write_gz, (saveto, lines))))
    return handles


scaling = normalization.median_of_ratios()

# Load effective genome size
with open(ld.EFFECTIVE_GENOME_SIZE, 'rb') as stream:
    effgsize = pickle.load(stream)
effgsize = {k: sum(v.values()) for k, v in effgsize.items()}

# Load all RNA models
with open(ld.RNA_MODELS, 'rb') as stream:
    rna_models = pickle.load(stream)

# Load all comparisons presets
presets = ld.Config.load()

handles = []
engine = rp.Reaper(threads=-1)
for cmp in presets:
    # Calculate scaling factors
    assembly = utils.assembly.get(organism=cmp.host)
    reads, _, scfactors = cmp.scaling(scaling[assembly.name])
    scfactors = {k: 1 / v for k, v in scfactors.items()}

    # Baseline
    seqlens = utils.assembly.seqsizes(cmp.organism)
    length = sum(seqlens.values())
    baseline = reads["control"].sum() / length

    # Construct the config
    print(f"[{cmp.ind}] Scaling factors: {scfactors}")
    model = cmp.model.set_control_baseline(baseline)
    enrichment = rp.cmp.Enrichment().set_scaling(scfactors["signal"], scfactors["control"])
    pcalling = cmp.pcalling
    nms = cmp.nms

    # Construct the workload
    workload = rp.Workload()
    for seq, length in seqlens.items():
        seqmodel, seqnms = copy.deepcopy(model), copy.deepcopy(nms)
        for strand in Strand.Forward, Strand.Reverse:
            rnas = rna_models[assembly.name].get((seq, strand), [])
            orientation = strand.to_orientation()
            if rnas:
                seqmodel.add_control_model(orientation, rnas, False, [256, 512, 1024])
                seqnms.add_regions(orientation, False, rnas)
            else:
                print(f"[{cmp.ind}] No RNA models for {seq} {strand}")

        config = rp.Config(seqmodel, enrichment, pcalling, seqnms)
        workload.add_region(seq, 0, length, config)

    # Add the comparison to the engine
    for tag, exps in (("Signal", cmp.signal), ("Control", cmp.control)):
        for exp in exps:
            source, layout = nfcore.rnaseq.extract.bam(exp, factory=reader)
            engine.add_source(tag, source, layout)

    engine.add_comparison(cmp.ind, "Signal", "Control", workload)
    harvest = engine.run()
    assert len(harvest) == 1, (cmp.ind, harvest)
    harvest = harvest.pop()

    handles.extend(save_results(cmp, harvest.regions))

for ind, h in handles:
    h.get()
    assert h.successful(), f"Failed to save results for {ind} {h.get()}"
SAVERS.close()
