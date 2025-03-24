#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
  library(optparse)
})

# Define command-line options
option_list <- list(
  make_option(c("-c", "--counts"), type = "character", help = "Raw counts file from featureCounts", metavar = "file"),
  make_option(c("-s", "--samples"), type = "character", help = "Sample metadata file", metavar = "file"),
  make_option(c("-t", "--tpm"), type = "character", help = "TPM values file (genes × samples)", metavar = "file"),
  make_option(c("-o", "--output"), type = "character", help = "Output file path", metavar = "file"),
  make_option(c("-g", "--group"), type = "character", default = "group", help = "Group column name in metadata [default = group]"),
  make_option(c("--baseline"), type = "character", help = "Baseline group name"),
  make_option(c("--treatment"), type = "character", help = "Treatment group name")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Check required args
required <- c("counts", "samples", "tpm", "output", "treatment", "baseline")
missing <- required[!required %in% names(opt) | sapply(opt[required], is.null)]
if (length(missing) > 0) {
  stop(paste("Missing required arguments:", paste(missing, collapse = ", ")))
}


suppressPackageStartupMessages({
  library(DESeq2)
  library(tibble)
  library(dplyr)
})

# Load count matrix
counts <- read.delim(opt$counts, comment.char = "#")
counts <- counts[!grepl("^__", counts$Geneid), ]  # remove summary rows
rownames(counts) <- counts$Geneid
counts_mat <- counts[, grep("Geneid|Length", colnames(counts), invert = TRUE)]

# Load sample metadata
coldata <- read.delim(opt$samples, row.names = 1)
coldata[[opt$group]] <- as.factor(coldata[[opt$group]])

# Subset to relevant groups only
groups_to_use <- c(opt$baseline, opt$treatment)
samples_to_use <- rownames(coldata)[coldata[[opt$group]] %in% groups_to_use]
coldata <- coldata[samples_to_use, , drop = FALSE]
counts_mat <- counts_mat[, samples_to_use]

# Load TPM matrix
tpm <- read.delim(opt$tpm, row.names = 1)
tpm <- tpm[, samples_to_use]

# TPM mean calculation
meanTPM_baseline <- rowMeans(tpm[, coldata[[opt$group]] == opt$baseline, drop = FALSE])
meanTPM_treatment <- rowMeans(tpm[, coldata[[opt$group]] == opt$treatment, drop = FALSE])


# Calculate log2FC from TPM
epsilon <- 1e-6
log2FC <- log2((meanTPM_treatment + epsilon) / (meanTPM_baseline + epsilon))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(counts_mat),
                              colData = coldata,
                              design = as.formula(paste("~", opt$group)))
dds <- DESeq(dds)
res <- results(dds, contrast = c(opt$group, opt$treatment, opt$baseline))
padj <- res$padj
names(padj) <- rownames(res)

# Build output table
out <- data.frame(
  GeneID = rownames(tpm),
  MeanTPM_Baseline = meanTPM_baseline[rownames(tpm)],
  MeanTPM_Treatment = meanTPM_treatment[rownames(tpm)],
  FDR = padj[rownames(tpm)],
  Log2FoldChange = log2FC[rownames(tpm)]
)

# Write output
write.table(out, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE);

