```shell
pixi run setup/download-annotations
# pixi run setup/liftoff # Note: not required, as lifted annotations are included in the repository
pixi run setup/

```











This repository accompanies the publication *"Host cell Z-RNAs activate ZBP1 during virus infections"* and contains the
code required to reproduce the bioinformatics results presented in the paper.

The repository is organized into biological or technical "stories" and includes several setup scripts for downloading,
preprocessing, and indexing reference genomes. Each setup and story has a corresponding make target, allowing easy
execution via `make <target>` from the repository's root directory.

---

### Prerequisites

- [Rust](https://www.rust-lang.org/tools/install) (tested with v1.84.0)
- [CMake](https://cmake.org/) (tested with v3.31.6)
- [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) (tested with v2.0.2)
- HPC cluster with the [Slurm](https://slurm.schedmd.com/quickstart.html) scheduler

**Notes:**

* Rust is not included in the Micromamba environment to avoid dependency conflicts.
* System-wide installations or HPC modules for Micromamba and Rust have not been tested. Please install them locally for
  the specific user to avoid potential conflicts.
* Running Micromamba from inside a conda environment does not work. Ensure that all conda environments are deactivated
  before calling `make <target>`.

---

### Setup

1. Clone the repository and navigate to its root directory.
2. Run the following commands for the initial setup:

   ```bash
   # Prepare the Micromamba environments and install the REAT A-I editing analysis tool
   make setup/env
   
   # Download the reference genome and index the annotation.
   # Annotation (GFF) downloading is optional, as it is already included in the repository.
   make setup/annotation
   
   # Derive the pre-mapping library.
   # This step is optional; results are already included in the repository.
   # make setup/pre-mapping 
   ```

   These commands can be executed on either a login node or a compute node in a single-core setup. However, subsequent
   steps requiring genome indexing must be performed on a compute node with multiple cores and additional memory.

3. Request a compute node and perform genome indexing:

   ```bash
   # Request a compute node
   srun --account $(whoami) --job-name nfcore-indices --cpus-per-task 32 --mem-per-cpu 8G --pty bash -i

   # Run genome indexing
   make setup/nfcore-indices
   ```

---

### Nextflow

#### Data Download

Public sequencing data can be downloaded using the [`nf-core/fetchngs`](https://nf-co.re/fetchngs) pipeline. For
example, to download data for *PRJNA256013*:

```bash
micromamba activate zdott
cd stories/nfcore/rnaseq/SRA/PRJNA256013/fastq
nfcore run \
    -r 7544cb9297a0db754120bd1cb8d7df4586a60610 \
    nf-core/fetchngs \
    --input id.csv \
    --outdir results \
    -profile singularity

# Move FASTQ files to the target directory
mv results/fastq/*.fastq.gz .

# Optionally, verify MD5 checksums. Note that committed checksums may not match downloaded ones due to 
# data storage and handling variations in SRA/ENA.
# cat MD5.txt | xargs -I{} -P $(nproc) sh -c 'echo "{}" | md5sum -c -'

# Clean up unnecessary files
rm -rf work results .nfcore*
```

Similarly, data can be downloaded for other public sequencing projects in the `stories/nextflow/series/SRA folder`.
Newly generated data in `stories/nextflow/series/internal` will be released upon publication and made equivalently
accessible via SRA/ENA.

---

#### Data Processing

A custom fork of the [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq) pipeline processes raw sequencing data,
including alignment, quality control, and RNA abundance estimation. This fork features minor modifications, such as
additional pre-mapping functionality and improved resource management.

To process each sequencing experiment, navigate to its corresponding directory under `stories/nextflow/series/SRA` or
`stories/nextflow/series/internal` and run:

```bash
micromamba activate zdott
nfcore run \
    -c resources/nfcore.config \
    -params-file resources/params.yaml \
    -resume \
    -profile slurm,latency \
    -r 41e95d6a4a24d1f4fe4f2c50bb4d9a4744158c9b \
    alnfedorov/rnaseq

# Optionally, clean up the directory
# rm -rf work .nfcore*
```

All experiments in `series/SRA` and `series/internal` directories should be processed before proceeding with the
downstream analyses.

---

### Stories

Each story can be executed by running `make <story>` from the root directory. As with genome indexing, these stories
should be run on a compute node with multiple cores (~16) and ~16GB of RAM per core. However, most analyses require
significantly fewer resources.

```bash
make stories/annotation     # Annotation filtering and indexing
make stories/normalization  # Genome binning and read counting for normalization/QC [running]
make stories/qc             # Additional quality control analyses

make stories/HSV1           # HSV1-specific analysis
make stories/A2I            # A-to-I editing analysis

# RIP analysis
make stories/RIP/pcalling   # Peak calling
make stories/RIP/clustering # dsRNA prediction and clustering (requires 96 cores and 16GB RAM per core)
make stories/RIP/annotation # Annotation of predicted dsRNA clusters
make stories/RIP/plots      # RIP analysis visualization

# Aberrantome analysis
make stories/aberrantome/calculate # Statistical tests for aberrant transcription events
make stories/aberrantome/plot      # Visualization of aberrantome analysis
```

The results for each story will be saved in the `<story>/ld/results` directory. The `ld` directory stands for "local
data" and is used for storing intermediate files and final outputs.

**Note:** Minimum Free Energy (MFE) plots for HSV-1 RNAs are not generated by default, as they require additional
dependencies and TEX backend configuration. Refer to `stories/HSV1/plot-mfe.py` for details.
