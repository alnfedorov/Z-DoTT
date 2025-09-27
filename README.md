This repository accompanies the publication *"Host cell Z-RNAs activate ZBP1 during virus infections"* and contains the
code required to reproduce the bioinformatics results presented in the paper.

The current version of the repository includes a number of quality of life improvements and simplified setup
compared to the originally published codebase
available [here](https://github.com/alnfedorov/Z-DoTT/tree/83b3746ef787226e2cc34804619854c5105979f4).
Because of that, some minor differences in the results should be expected, which, however, do not affect the overall
conclusions of the study.

The repository is organized into five broad sections:

* `setup/`: Scripts for preparing essential resources, e.g., downloading annotations and indexing reference genomes.
* `assemblies/`: Reference genome assemblies, including fasta files, annotations, and Python annotation indexes.
* `data/`: Sequencing data used in the study and its metadata, including both public datasets and newly generated
  results.
* `utils/`: Utility scripts for various tasks used throughout the repository.
* `analyses/`: Main analysis scripts and workflows organized into separate, mostly independent modules.

All dependencies and execution targets are managed via [Pixi](https://pixi.sh/latest/) with all relevant tasks wrapped
into a separate pixi command (see below).

---

### Prerequisites

- [Pixi](https://pixi.sh/latest/) binary (tested on v0.55.0).
- An x86-64 HPC cluster with the [Slurm](https://slurm.schedmd.com/quickstart.html) scheduler running a
  Linux-based OS (tested on Ubuntu 22.04.5).

**Notes:**

* To avoid environment conflicts, ensure that all conda/mamba environments are deactivated before running any pixi
  commands.
* All pixi tasks are *naive* and thus will re-run the entire task even if no code or data has changed. This is
  intentional as setting up proper caching would be overly complex for this repository.

---

### Setup

1. Clone the repository and navigate to its root directory:

    ```shell
    git clone https://github.com/alnfedorov/Z-DoTT.git
    cd Z-DoTT
    ```

2. Run the following commands for the initial setup:

    ```shell
    # Install all dependencies via Pixi
    pixi install --all --frozen
    
    # Download genome assemblies and make Python indexes
    pixi run setup/download-annotations
    pixi run setup/index-annotations
   
    # Compile REAT, a tool for RNA editing analysis
    pixi run setup/REAT
    
    # Note: These steps are optional, as the required annotations are already included in the repository.
    # pixi run setup/liftoff # lift over annotations from GRCh38 to CHM13v2 (min 64GB RAM required)
    # pixi run setup/derive-premap-library # derive the pre-mapping rRNA library
    ```

   These commands can be executed on either a login or a compute node as they do not require significant resources,
   except for the optional liftoff.

3. Request a compute node and make STAR/Salmon indexes for (pseudo)alignment of sequencing data:

    ```shell
    srun --account $(whoami) --job-name nfcore-indices \
    --cpus-per-task 24 --mem-per-cpu 8G --time 08:00:00 --pty \
    bash -lc "pixi run setup/make-nfcore-indexes"
    ```

---

### Sequencing data

#### Download

Sequencing data can be downloaded using the [`nf-core/fetchngs`](https://nf-co.re/fetchngs) pipeline. For
example, to download data for *PRJNA256013*:

```shell
# Run on any network-connected node (login or compute)
pixi run nfcore/fetchngs PRJNA256013

# List downloaded files
ls -alh data/PRJNA256013/fastq
```

Similarly, raw data can be downloaded for other sequencing projects listed in the `data` folder (e.g., B256178). In rare
instances, `nf-core/fetchngs` may fail due to network issues. A simple re-run usually resolves the problem.

#### Processing and QC

A custom fork of the [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq) pipeline is used to process raw sequencing
data, including alignment, quality control, and RNA abundance estimation. This fork features minor modifications, such
as additional pre-mapping functionality and improved resource management.

Like `fetchngs`, the RNA-seq pipeline can be executed via a single pixi command. For example, to process the data for
*PRJNA256013*:

```shell
# Run on a *login* node
pixi run nfcore/rnaseq PRJNA256013

# List processed files
ls -alh analyses/nfcore/rnaseq/PRJNA256013/results

# Optionally, clean up the caches
# rm -rf analyses/nfcore/rnaseq/PRJNA256013/{work,.nfcore*}
```

The pipeline should be run on a login node, as it submits jobs to the cluster via the Slurm scheduler. If needed,
users can tweak the parameters using the `params.yaml` and `nextflow.config` files located in the
`analyses/nfcore/resources` directory. Refer to the [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq) documentation
for details.

All experiments in the `analyses/nfcore/rnaseq` directory should be processed before proceeding with the
downstream analyses.

---

### Analyses

Most analyses feature a dedicated `pixi` task that will execute the entire (sub)analysis from start to finish. Usually,
each analysis is organized as follows:

```text
analyses/
└── example                 # Example analysis folder
    │   __init__.py         # Each folder is a Python package and optionally contains an __init__.py file
    ├── resources           # General resources required for the analysis, including local utils and data files
    │   ├── __init__.py
    │   ├── utils.py
    │   └── regions.bed
    ├── results             # Final or intermediate results, including caches, are stored here
    │   ├── __init__.py
    │   ├── cache.pkl
    │   └── result-1.csv
    ├── script-1.py         # Main analysis script, almost exclusively in Python
    └── script-2.py
```

Note, that most `resources`, `results`, and most analysis folders are Python packages (i.e., they contain an
`__init__.py` file). These files usually export paths to relevant subdirectories, but `resources` and analysis package
also re-export functions, data models, or configs (particularly in `resources`).

In rare cases, analyses are split into subanalyses, e.g., `stories/RIP` and `stories/aberrantome`. In such cases,
each subanalysis is still organized as above and the grouping is largely for convenience.

As with STAR/Salmon indexing, analyses should be run on a compute node with multiple cores (~16) and ~8GB of RAM
per core. However, most analyses require significantly fewer resources. Notable exceptions are detailed below.

Full list of available analyses:

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

You can find plots presented in the manuscript in the respective `analyses/<name>/results` subfolders.

**Note:** Minimum Free Energy (MFE) plots for HSV-1 RNAs are not generated by default, as they require additional
dependencies and TEX backend configuration. Refer to `analyses/HSV1/plot-mfe.py` for details.
