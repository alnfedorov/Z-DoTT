from pathlib import Path

from joblib import Parallel, delayed
from pybedtools import BedTool

import ld
import utils
from assemblies import GRCm39, CHM13v2
from stories import annotation


def job(assembly):
    # 1. Annotate target editing sites
    anno = annotation.load.gencode(assembly.name)
    categories, mapping = ld.utils.annotate(assembly, anno, ld.c.tracks.all_passed / f"{assembly.name}.bed.gz")

    # Save the mapping
    saveto = ld.c.tracks.all_annotated / assembly.name
    saveto.mkdir(parents=True, exist_ok=True)

    for location, sites in categories.items():
        utils.bed.tbindex(BedTool(sites).sort(), saveto / f"{location}.bed.gz")

    # 2. Annotate each sample
    def process(dataset: str, sample: Path):
        df = ld.utils.load.filtered(sample)

        category = []
        for contig, strand, pos in zip(df['contig'], df['trstrand'], df['pos']):
            category.append(mapping[contig, strand, pos])
        df['location'] = category

        saveto = ld.c.reat.annotated / assembly.name / dataset
        saveto.mkdir(parents=True, exist_ok=True)

        df.to_csv(saveto / sample.name, index=False)

    workload = []
    for dataset in (ld.c.reat.filtered / assembly.name).iterdir():
        for sample in dataset.iterdir():
            workload.append((dataset.name, sample))

    Parallel(n_jobs=-1)(delayed(process)(*args) for args in workload)


Parallel(n_jobs=-1, backend='threading')(delayed(job)(assembly) for assembly in [GRCm39, CHM13v2])
