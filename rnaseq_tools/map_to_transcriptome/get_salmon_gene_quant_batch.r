#!/usr/bin/env Rscript

# Get arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript get_salmon_gene_quant_batch.r <sample_dirs.txt> <tx2gene.tsv> <out_prefix>")
}

# Load libraries
suppressMessages(library(tximport))
suppressMessages(library(readr))

sample_list_file <- args[1]
tx2gene_file <- args[2]
out_prefix <- args[3]

# Read list of sample directories
sample_dirs <- readLines(sample_list_file)

# Create named vector of quant.sf file paths
quant_files <- file.path(sample_dirs, "quant.sf")
names(quant_files) <- basename(sample_dirs)

# Check all quant.sf files exist
missing <- !file.exists(quant_files)
if (any(missing)) {
  stop("Missing quant.sf files in the following directories:\n",
       paste(sample_dirs[missing], collapse = "\n"))
}

# Read transcript-to-gene mapping
tx2gene <- read_tsv(tx2gene_file, col_names = TRUE)

# Run tximport
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")

# Save counts, TPMs, and lengths
write.table(data.frame(GeneID=rownames(txi$counts), txi$counts), file= paste0(out_prefix, ".gene_counts.tsv"), sep= "\t", quote= FALSE, row.names= FALSE)
write.table(data.frame(GeneID=rownames(txi$abundance), txi$abundance), file = paste0(out_prefix, ".gene_tpm.tsv"), sep= "\t", quote= FALSE, row.names= FALSE)
write.table(data.frame(GeneID=rownames(txi$length), txi$length), file = paste0(out_prefix, ".gene_lengths.csv"), sep= "\t", quote= FALSE, row.names= FALSE)
saveRDS(txi, paste0(out_prefix, ".RDS")); # This is kept for DESeq2 analysis.

