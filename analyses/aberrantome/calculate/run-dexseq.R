library(DEXSeq)

# Parallelization
# library("BiocParallel")
# BPPARAM = MulticoreParam(multicoreWorkers())
# register(BPPARAM)

args <- commandArgs(TRUE)

SAMPLES <- args[1]
COUNTS <- args[2]
SAVETO <- args[3]
BASELINE <- args[4]
NORMALIZED <- as.logical(args[5])

# SAMPLES <- '/tmp/tmps54gewnj'
# COUNTS <- '/tmp/tmpkk5l5m4h'
# SAVETO <- '/tmp/tmpxfqyuv7d'
# BASELINE <- 'control'
# NORMALIZED <- as.logical('True')

samples <- read.csv(SAMPLES, header = TRUE, row.names = 1)
samples$condition <- relevel(factor(samples$condition), ref = BASELINE)

counts <- read.csv(COUNTS, header = TRUE)
groups <- counts$group
counts$group <- NULL

elements <- counts$element
counts$element <- NULL

stopifnot(rownames(samples) == colnames(counts))

dxd <- DEXSeqDataSet(
  countData = counts,
  sampleData = samples,
  design = ~sample + exon + condition:exon,
  featureID = elements,
  groupID = groups
)

if (NORMALIZED == FALSE) {
  dxd <- estimateSizeFactors(dxd)
} else {
  normFactors <- rep(1, nrow(samples))
  # I have no idea why this is necessary. Adopted from the DEXSeq source code.
  sizeFactors(dxd) <- rep(normFactors, 2)
  maxExons <- length(unique(dxd@modelFrameBM$exon))
  dxd@modelFrameBM$sizeFactor <- rep(1, each = maxExons)
}
print("DEXSeq size factors:")
print(sizeFactors(dxd))

dxd <- estimateDispersions(dxd, quiet = TRUE)
dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd)

dxr <- DEXSeqResults(dxd)

dxr$countData <- NULL
dxr$genomicData <- NULL

write.csv(dxr, gzfile(SAVETO))
