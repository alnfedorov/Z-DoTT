This repository accompanies the publication *"Host cell Z-RNAs activate ZBP1 during virus infections"* and contains the
code required to reproduce the bioinformatics results presented in the paper.

This version introduces several quality-of-life improvements and a simplified setup compared to the original codebase (
available [here](https://github.com/alnfedorov/Z-DoTT/tree/83b3746ef787226e2cc34804619854c5105979f4). As a result, you
may observe minor differences in outputs; these do not affect the study’s overall conclusions.

This repository is provided *as is* for reference purposes only. Please open a GitHub issue if you encounter bugs or
problems reproducing our results. However, I can’t assist with repurposing the code for other studies, as the current
methods are highly experimental and extensively tailored to our datasets.

-----

## Repository Structure

The repository is composed of five main folders:

* `setup`: Scripts for preparing essential resources, like downloading annotations and indexing reference genomes.
* `assemblies`: Reference genome assemblies, including FASTA files, annotations, and Python annotation indexes.
* `data`: Raw sequencing data and metadata used in the study.
* `utils`: Utility scripts used throughout the repository.
* `analyses`: The main analysis scripts and workflows, organized into distinct modules.

The entire workflow is managed using [Pixi](https://pixi.sh/latest/); commands are exposed as `pixi run <task>`.

---

## Prerequisites

* [Pixi](https://pixi.sh/latest/) (tested with v0.55.0).
* Access to an x86_64 HPC cluster running Linux (tested on Ubuntu 22.04.5) with
  the [Slurm](https://slurm.schedmd.com/quickstart.html) scheduler.

**Important Notes:**

* To prevent environment conflicts, please ensure any active conda or mamba environments are **deactivated** before
  running Pixi commands.
* Pixi tasks are cacheless and re-execute from scratch each time they are run.
* Pay attention to exit codes and error messages. If a task fails, it will not produce the expected output, and
  subsequent tasks depending on it will also fail.

---

### Setup

1. Clone the repository and navigate into the root directory:

    ```shell
    git clone https://github.com/alnfedorov/Z-DoTT.git
    cd Z-DoTT
    ```

2. Perform the initial setup:

    ```shell
    # Install all dependencies from the frozen lockfile
    pixi install --all --frozen
    
    # Download genome assemblies and create Python indexes
    pixi run setup/download-annotations
    pixi run setup/index-annotations
   
    # Compile REAT, a tool for RNA editing analysis
    pixi run setup/REAT
    
    # The following steps are optional as their outputs are already included in the repo
    # pixi run setup/liftoff # Lifts over annotations from GRCh38 to CHM13v2 (requires ~80GB RAM)
   srun --partition long,short --account $(whoami) --job-name biobit 
   --cpus-per-task 8 --mem-per-cpu 10G --time 0-08:00:00 --pty bash -i
    # pixi run setup/derive-premap-library # Derives the pre-mapping rRNA library
    ```

   These commands can be run on a login node, as they are not resource-intensive (with the exception of the optional
   `liftoff` step).

3. Request a compute node to generate the STAR and Salmon indexes required for sequence alignment:

    ```shell
    srun --account $(whoami) --job-name nfcore-indices \
    --cpus-per-task 18 --mem-per-cpu 8G --time 04:00:00 --pty \
    bash -lc "pixi run setup/make-nfcore-indexes"
    ```

---

### Sequencing data

#### Download

Sequencing data can be downloaded from public archives using the [`nf-core/fetchngs`](https://nf-co.re/fetchngs)
pipeline.

For example, to download data for *PRJNA256013*:

```shell
# Run on any network-connected node
pixi run nfcore/fetchngs PRJNA256013

# List downloaded files
ls -alh data/PRJNA256013/fastq
```

This should be repeated for all other projects listed in the `data/` directory. If a download fails due to network
issues, re-running the command usually resolves the problem.

#### Processing and QC

A custom fork of the [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq) pipeline is used to process raw sequencing
data. This fork includes minor modifications for pre-mapping and resource management.

To process the *PRJNA256013* dataset, run the following command from a **login node**. The pipeline will automatically
submit jobs to the Slurm cluster.

```shell
# Run on a *login* node
pixi run nfcore/rnaseq PRJNA256013

# List processed files
ls -alh workflows/nfcore/rnaseq/PRJNA256013/results

# Optionally, clean up the caches
# rm -rf workflows/nfcore/rnaseq/PRJNA256013/{work,.nfcore*}
```

Pipeline parameters can be modified in the `params.yaml` and `nextflow.config` files located in
`analyses/nfcore/resources`. For more details, refer to the official
`[nf-core/rnaseq documentation](https://github.com/nf-core/rnaseq\)`.

**Important:** All datasets in `analyses/nfcore/rnaseq` must be downloaded and processed before running the downstream
analyses.


---

### Analyses

Most analyses have a dedicated `pixi` task that executes the entire workflow. A typical analysis directory looks like:

```text
analyses/
└── example/                # Example analysis folder
    ├── resources/          # Resources like local utils and data files
    ├── results/            # Final and intermediate results, including caches
    └── script-1.py         # Main analysis script(s)
```

Note that `resources`, `results`, and most analysis folders are commonly organized as Python packages (i.e., they
contain an `__init__.py` file). These files export paths to relevant subdirectories; `resources` and the root analysis
package may also re-export functions, data models, or configs.

In rare cases, analyses are split into subanalyses, e.g., `stories/RIP` and `stories/aberrantome`. Each sub-analysis
retains the structure above; the grouping is for convenience.

Most analyses should be run on a multi-core compute node (e.g., 16 cores with \~8GB RAM per core), although many require
fewer resources. Notable exceptions with higher requirements are listed below.

Order of execution:

```shell
# Annotation filtering and indexing
pixi run workflows/annotation

# Genome binning and read counting for QC
pixi run workflows/normalization

# Additional quality control workflows
pixi run workflows/qc

# HSV-1-specific analysis
pixi run workflows/HSV1

# A-to-I editing analysis
pixi run workflows/A2I

# ---- RIP Analysis ----
# Peak calling
pixi run workflows/RIP/pcalling
# dsRNA prediction and clustering (requires 96 cores & 16GB RAM/core)
pixi run workflows/RIP/clustering
# Annotation of predicted dsRNA clusters
pixi run workflows/RIP/annotation
# Visualization of RIP analysis results
pixi run workflows/RIP/plots

# ---- Aberrantome Analysis ----
# Statistical tests for aberrant transcription events
pixi run workflows/aberrantome/calculate
# Visualization of aberrantome analysis results
pixi run workflows/aberrantome/plot
```

The final plots presented in the manuscript can be found in the respective `analyses/<name>/results` subfolders.

**Note:** Minimum Free Energy (MFE) plots for HSV-1 RNAs are not generated by default, as they require additional
dependencies and TEX backend configuration. Refer to `analyses/HSV1/plot-mfe.py` for details.
