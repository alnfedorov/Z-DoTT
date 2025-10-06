library("DESeq2")

# Parallelization
library("BiocParallel")
register(MulticoreParam(multicoreWorkers()))

args = commandArgs(TRUE)

SAMPLES = args[1]
COUNTS = args[2]
DESIGN = args[3]
ATTRIBUTE = args[4]
TARGET = args[5]
REFERENCE = args[6]
ALTERNATIVE = args[7]
LOG2FC = args[8]
PADJ = args[9]
SAVETO = args[10]

LOG2FC = as.numeric(LOG2FC)
PADJ = as.numeric(PADJ)
DESIGN = as.formula(DESIGN)

# Prepare samples table
samples = read.csv(SAMPLES, header = TRUE, row.names = 1)
samples[, ATTRIBUTE] = relevel(factor(samples[, ATTRIBUTE]), ref = REFERENCE)

# Load counts
counts = read.csv(COUNTS, header = TRUE, row.names = 1)

# All samples must be present in the counts table
stopifnot(rownames(samples) == colnames(counts))

# Make the test
dds = DESeqDataSetFromMatrix(countData = counts, colData = samples, design = DESIGN)
sizeFactors(dds) = rep(1, ncol(dds))
dds = DESeq(dds)

name = paste(ATTRIBUTE, paste(TARGET, REFERENCE, sep="_vs_"), sep = "_")
result = results(dds, name = name, alpha = PADJ, lfcThreshold = LOG2FC, altHypothesis = ALTERNATIVE)
result = lfcShrink(dds, coef = name, res = result, type = "apeglm")

write.csv(result, gzfile(SAVETO))
