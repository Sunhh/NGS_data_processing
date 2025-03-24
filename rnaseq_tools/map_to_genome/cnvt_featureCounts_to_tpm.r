#!/usr/bin/env Rscript

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript calculate_tpm.R <featureCounts_file> <output_file>")
}

input_file <- args[1]
output_file <- args[2]

# Read featureCounts output
fc <- read.delim(input_file, comment.char = "#")

# Remove summary rows if present (e.g., "__no_feature")
fc <- fc[!grepl("^__", fc[,1]), ]

# Extract gene lengths and raw counts
gene_lengths <- fc$Length  # column usually named 'Length'
raw_counts <- fc[, 7:ncol(fc)]  # columns after column 6 are sample counts

# Calculate RPK
rpk <- raw_counts / (gene_lengths / 1000)

# Calculate TPM
tpm <- apply(rpk, 2, function(x) x / sum(x) * 1e6)

# Combine with gene ID
tpm_df <- data.frame(GeneID = fc$Geneid, tpm)

# Write to output
write.table(tpm_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
