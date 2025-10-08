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
