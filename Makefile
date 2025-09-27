# Run all commands in a single shell
.ONESHELL: *
# All targets are phony
.PHONY: $(MAKECMDGOALS)

# Use bash syntax everywhere
SHELL := $(shell which bash)
# Interactive sessions (to properly register micromamba) with extra strictness (fail on errors/pipefails)
.SHELLFLAGS := -i -e -o pipefail -c
# Delete target if command fails
.DELETE_ON_ERROR:

# Warn if an undefined variable is referenced
MAKEFLAGS += --warn-undefined-variables
# Disable default suffix rules
MAKEFLAGS += --no-builtin-rules

# Common setup command for activating the environment
define ACTIVATE_ENV
micromamba activate zdott
export PYTHONPATH="$$(pwd)"
endef

setup/pre-mapping:
	$(ACTIVATE_ENV)
	cd setup/pre-mapping
	python derive-premap-library.py

setup/nextflow-indices:
	$(ACTIVATE_ENV)
	ASSEMBLIES="$$(pwd)/assemblies"
	PREMAPPING="$$(pwd)/setup/pre-mapping/results"
	cd stories/nextflow/resources

	mkdir -p indexes/GRCm39 indexes/CHM13v2

	# Cat the annotation data
	zcat "$${ASSEMBLIES}/GRCm39/GRCm39.primary_assembly.genome.fa.gz" > indexes/GRCm39/sequence.fa
	zcat "$${ASSEMBLIES}/GRCm39/gencode/gencode.vM36.primary_assembly.annotation.gff3.gz" > indexes/GRCm39/annotation.gff

	zcat "$${ASSEMBLIES}/CHM13v2/CHM13v2.fa.gz" > indexes/CHM13v2/sequence.fa
	zcat "$${ASSEMBLIES}/CHM13v2/gencode/CHM13v2.liftoff+gencode-47.gff3.gz" > indexes/CHM13v2/annotation.gff

	for assembly in GRCm39 CHM13v2;
	do
		# Add viruses
		for virus in HSV1 IAV;
		do
			cat "$${ASSEMBLIES}/$${virus}/sequence.fasta" >> "indexes/$${assembly}/sequence.fa"
			cat "$${ASSEMBLIES}/$${virus}/sequence.gff3" >> "indexes/$${assembly}/annotation.gff"
		done

		# Generate pre-mapping index
		bash premap-index.sh "indexes/$${assembly}" "$${PREMAPPING}/$${assembly}"

		# Run the main indexing
		bash main-index.sh "indexes/$${assembly}"

		# Link configs
		for fname in nextflow.config params.yaml;
		do
			ln -s "$$(pwd)/$${fname}" "indexes/$${assembly}/$${fname}"
		done
	done

stories/annotation:
	$(ACTIVATE_ENV)
	cd stories/annotation
	python filter-gencode.py
	python index-rna-boundaries.py && python derive-rna-cores.py &
	python resolve-gencode.py &
	wait

stories/normalization:
	$(ACTIVATE_ENV)
	cd stories/normalization
	python bin-genomes.py
	python count-reads.py
	python vRNA-load.py

stories/qc:
	$(ACTIVATE_ENV)
	cd stories/qc
	python count-biotypes.py
	python plot-pca.py
	python plot-biotypes.py

stories/HSV1:
	$(ACTIVATE_ENV)
	cd stories/HSV1

	# Make sure that we have the required comparisons info
	python ../RIP/pcalling/make-config.py

	# Experimental data
	python calculate-enrichment.py
	python prepare-annotations.py

	python plot-circos.py
	python plot-enrichment-ratio.py
	python plot-vRNA-ratios.py

	# Min-free-energy predictions
	python predict-mfe.py
	python zh-score-mfe.py
	# python plot-mfe.py # Disabled because it requires extra packages. Run manually if needed

stories/A2I:
	$(ACTIVATE_ENV)
	cd stories/A2I

	python candidates-generation.py
	python candidates-filtering.py
	python candidates-annotation.py

	python plot-time-series.py &
	python plot-bars.py &
	wait

stories/RIP/pcalling:
	$(ACTIVATE_ENV)
	cd stories/RIP/pcalling
	python make-config.py
	python make-rna-models.py
	python call-peaks.py

stories/RIP/clustering:
	$(ACTIVATE_ENV)
	cd stories/RIP/clustering

	python make-config.py
	python peaks-prefiltering.py
	python derive-insulators.py

	# dsRNA prediction and filtering
	python dsRNA-prediction.py
	python dsRNA-filtering.py

	# Final peaks filtering and co-clustering with passed dsRNAs
	python peaks-filtering.py
	python clustering.py

stories/RIP/annotation:
	$(ACTIVATE_ENV)
	cd stories/RIP/annotation

	python make-config.py

	# Cache counts and the signal
	python cache-signal.py &
	python cache-counts.py &
	wait

	# Annotate all peaks/dsRNAs
	python A2I-editing.py &
	python dsRNA-gaps.py &
	python qPCR-probes.py &
	python localization.py &
	python sequence-composition.py &
	python stat-comparisons.py &
	wait

	# Make final summaries
	python make-summaries.py

stories/RIP/plots:
	$(ACTIVATE_ENV)
	cd stories/RIP/plots

	python resolve.py
	python venn-diagram.py &
	python seqloc-distribution.py &
	python loop-size.py &
	python volcano-plot.py &
	python editing-summary.py &
	wait

stories/aberrantome/calculate:
	$(ACTIVATE_ENV)
	cd stories/aberrantome/calculate

	python make-config.py
	python prepare-rnas.py
	python count-reads.py
	python run-dexseq.py
	python annotate-results.py

stories/aberrantome/plot:
	$(ACTIVATE_ENV)
	cd stories/aberrantome/plot

	python resolve.py
	python scores-regplot.py &
	python summary.py &
	python volcano-plot.py &
	python z-rna-concordance.py &
	wait
