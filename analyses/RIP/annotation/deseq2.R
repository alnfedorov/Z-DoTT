library("DESeq2")

# Parallelization
library("BiocParallel")
register(MulticoreParam(multicoreWorkers()))

args = commandArgs(TRUE)
DESIGN = args[1]
COUNTS = args[2]
SIZE_FACTORS = args[3]
TEST_FORMULA = args[4]
BASELINES = unlist(strsplit(args[5], "$", fixed = TRUE))
COMPARISONS = unlist(strsplit(args[6], "$", fixed = TRUE))
ALTERNATIVE = args[7]
LOG2FC = as.numeric(args[8])
ALPHA = as.numeric(args[9])
SAVETO = args[10]

# Load the design table
colData = read.table(DESIGN, header=TRUE, row.names=1, sep="\t")
for (bs in BASELINES) {
    bs = unlist(strsplit(bs, "%", fixed = TRUE))
    colData[, bs[1]] = relevel(factor(colData[, bs[1]]), ref = bs[2])
}

# Load the counts
counts = read.table(COUNTS, header = TRUE, row.names = "index", sep = "\t")
stopifnot(rownames(colData) == colnames(counts))

# Make the formula
formula = as.formula(TEST_FORMULA)

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = colData,
  design = formula
)

# Remove low count records
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep,]

# Size factors
factors = read.table(SIZE_FACTORS, header = TRUE, sep = "\t", row.names = 1)
stopifnot(rownames(samples) == colnames(factors))
sizeFactors(dds) = as.numeric(as.vector(factors[1,]))

# Run the test
dds = DESeq(dds)
print(resultsNames(dds))

for (test in COMPARISONS) {
    test = unlist(strsplit(test, "%", fixed = TRUE))
    saveto = test[1]
    name = test[2]

    res = results(dds, name = name, alpha = ALPHA, lfcThreshold = LOG2FC, altHypothesis = ALTERNATIVE)
    res = lfcShrink(dds, coef=name, type="apeglm", res=res)

    # Save as .csv.gz files
    saveto = file.path(SAVETO, paste(saveto, "csv.gz", sep="."))
    write.csv(res, gzfile(saveto))
}
