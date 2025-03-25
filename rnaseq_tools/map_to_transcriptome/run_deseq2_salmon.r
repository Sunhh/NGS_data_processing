#!/usr/bin/env Rscript
# Ref: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#DESeq2
# Ref: https://combine-lab.github.io/salmon/getting_started/
# Ref: https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

suppressMessages(library(optparse))
# Define command line options
option_list <- list(
  make_option(c("-r", "--rds"), type = "character", help = "Input RDS file from tximport", metavar = "FILE"),
  make_option(c("-m", "--metadata"), type = "character", help = "Sample metadata CSV file with columns: sample, group", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", help = "Output file name for DE results", metavar = "FILE"),
  make_option(c("-b", "--baseline"), type = "character", help = "Baseline group name (e.g., control)"),
  make_option(c("-t", "--treatment"), type = "character", help = "Treatment group name (e.g., treated)")
)

suppressMessages(library(tximport))
suppressMessages(library(DESeq2))
suppressMessages(library(readr))
suppressMessages(library(dplyr))

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$rds) || is.null(opt$metadata) || is.null(opt$output) ||
    is.null(opt$baseline) || is.null(opt$treatment)) {
  print_help(opt_parser)
  stop("All arguments are required.\n", call. = FALSE)
}

# Load tximport object and sample metadata
txi <- readRDS(opt$rds)
samples <- read.table(opt$metadata, sep= "\t", header= TRUE, stringsAsFactors = FALSE)

# Check required columns
if (!all(c("sample", "group") %in% colnames(samples))) {
  stop("Metadata file must contain 'sample' and 'group' columns.")
}

## Ensure sample order matches
#samples <- samples[match(colnames(txi$counts), samples$sample), ]
#if (any(is.na(samples$sample))) {
#  stop("Mismatch between sample names in metadata and quantification files.")
#}
#rownames(samples) <- samples$sample

# Filter for baseline and treatment samples
samples <- samples[samples$group %in% c(opt$baseline, opt$treatment), ]

# Set rownames and match sample order
rownames(samples) <- samples$sample
samples <- samples[intersect(colnames(txi$counts), samples$sample), ]
if (nrow(samples) == 0) {
  stop("No samples belong to the specified baseline or treatment groups.")
}

# Subset txi object
keep_samples <- samples$sample
txi$counts <- txi$counts[, keep_samples, drop = FALSE]
txi$abundance <- txi$abundance[, keep_samples, drop = FALSE]
txi$length <- txi$length[, keep_samples, drop = FALSE]

# Convert group to factor and set baseline
samples$group <- factor(samples$group, levels = c(opt$baseline, opt$treatment))

# Safety checks
n_base <- sum(samples$group == opt$baseline)
n_treat <- sum(samples$group == opt$treatment)
if (n_base < 2 || n_treat < 2) {
  stop(paste("Need at least 2 samples in each group. Found:",
             n_base, opt$baseline, "and", n_treat, opt$treatment))
}

# Run DESeq2
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ group)
dds <- DESeq(dds)

# Run the comparison
res <- results(dds, contrast = c("group", opt$treatment, opt$baseline))
res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)

# Compute mean TPM per group
abund <- txi$abundance
grouped_samples <- split(colnames(abund), samples$group)
mean_tpm_baseline  <- rowMeans(abund[, grouped_samples[[opt$baseline]], drop = FALSE])
mean_tpm_treatment <- rowMeans(abund[, grouped_samples[[opt$treatment]], drop = FALSE])

# Add mean TPMs and rename columns
res_df$MeanTPM_Baseline  <- mean_tpm_baseline[res_df$GeneID]
res_df$MeanTPM_Treatment <- mean_tpm_treatment[res_df$GeneID]

output_df <- res_df %>%
  select(GeneID, MeanTPM_Baseline, MeanTPM_Treatment, 
         FDR = padj, Log2FoldChange = log2FoldChange)

# Save to TSV
write.table(output_df, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE)

