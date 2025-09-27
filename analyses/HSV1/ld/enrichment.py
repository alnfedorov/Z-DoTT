from collections.abc import Iterable
from pathlib import Path
from subprocess import check_call

import numpy as np
import numpy.typing as npt
import pyBigWig
from attrs import define
from biobit.core.loc import PerStrand


@define
class EnrichmentTracks:
    signal: PerStrand[dict[str, npt.NDArray[np.float32]]]
    control: PerStrand[dict[str, npt.NDArray[np.float32]]]
    enrichment: PerStrand[dict[str, npt.NDArray[np.float32]]]
    scaled_enrichment: PerStrand[dict[str, npt.NDArray[np.float32]]]


def calculate(
        regions: dict[str, int], ind: str, cache: Path,
        signal: list[Path | str], signal_scale: float,
        control: list[Path | str], control_scale: float,
        pseudo_cnt: float = 1, cpus: int = 1,
        exclude_flags: int = 2572, include_flags: int = 3, minmapq: int = 0
) -> EnrichmentTracks:
    # 1. Create a bed record with target regions
    regions_bed = cache / f"{ind}.regions.bed"
    with open(regions_bed, "w") as stream:
        for name, length in regions.items():
            stream.write(f"{name}\t0\t{length}\n")

    # 2. Subsamble and merge target BAM files
    cleaned_files = {}
    for prefix, files in ("signal", signal), ("control", control):
        if len(files) == 1:
            target = files[0]
        else:
            target = cache / f"{ind}.{prefix}.merged.bam"
            check_call([
                "samtools", "merge", "-@", str(cpus), "-L", regions_bed, "--write-index",
                "-o", target, *files
            ])

        saveto = cache / f"{ind}.{prefix}.bam"
        check_call([
            "samtools", "view", "-@", str(cpus), "-L", regions_bed, "--write-index", "-b",
            "--excl-flags", str(exclude_flags), "--incl-flags", str(include_flags), "--min-MQ", str(minmapq),
            "-o", saveto, target
        ])
        cleaned_files[prefix] = saveto

    # 3. Calculate the coverage using deeptools bamCoverage
    scfactors = {"signal": signal_scale, "control": control_scale}
    coverage = {}
    for name, bam in cleaned_files.items():
        for postfix, strand in ("fwd", "forward"), ("rev", "reverse"):
            saveto = cache / f"{ind}.{name}.{postfix}.bigWig"
            check_call([
                "bamCoverage", "-b", bam, "-o", saveto, "--filterRNAstrand", strand,
                "--outFileFormat", "bigwig", "--binSize", "1", "--numberOfProcessors", str(cpus),
            ])

            values = {}
            with pyBigWig.open(saveto.as_posix()) as file:
                for region, length in regions.items():
                    v = file.values(region, 0, length, numpy=True)
                    v[np.isnan(v)] = 0
                    v *= scfactors[name]
                    values[region] = v.astype(np.float32)

            coverage[name, strand] = values

    # 4. Calculate the enrichment
    enrichment = {}
    for strand in "forward", "reverse":
        signal, control = coverage["signal", strand], coverage["control", strand]
        enrichment[strand] = {
            region: (signal[region] + pseudo_cnt) / (control[region] + pseudo_cnt) for region in regions
        }
        to_bigwig(enrichment[strand], cache / f"{ind}.enrichment.{strand}.bigWig")

    # 5. Calculate scaled enrichment
    scaled_enrichment = {}
    for strand in "forward", "reverse":
        values = {}
        for region in regions:
            scale = coverage["signal", strand][region] / np.mean(coverage["signal", strand][region])
            values[region] = enrichment[strand][region] * scale
        scaled_enrichment[strand] = values
        to_bigwig(scaled_enrichment[strand], cache / f"{ind}.scaled_enrichment.{strand}.bigWig")

    result = EnrichmentTracks(
        signal=PerStrand(forward=coverage["signal", "forward"], reverse=coverage["signal", "reverse"]),
        control=PerStrand(forward=coverage["control", "forward"], reverse=coverage["control", "reverse"]),
        enrichment=PerStrand(forward=enrichment["forward"], reverse=enrichment["reverse"]),
        scaled_enrichment=PerStrand(forward=scaled_enrichment["forward"], reverse=scaled_enrichment["reverse"])
    )

    return result


def to_bigwig(values: dict[str, Iterable[float]], saveto: Path):
    values = {k: list(map(float, v)) for k, v in values.items()}
    sizes = sorted((k, len(v)) for k, v in values.items())

    with pyBigWig.open(saveto.as_posix(), "w") as file:
        file.addHeader(sizes)
        for name, _ in sizes:
            vals = values[name]

            seqid = [name] * len(vals)
            starts = list(range(0, len(vals)))
            ends = list(range(1, len(vals) + 1))
            file.addEntries(seqid, starts, ends=ends, values=vals)
