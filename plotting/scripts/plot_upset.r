#!/usr/bin/env Rscript
suppressMessages({
library(UpSetR)
library(tidyverse)
});

# Function to read gene lists from input files
read_gene_list <- function(file) {
  readLines(file) %>% unique()
}

# Main script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript script.R out_prefix <gene_list_1.txt> <gene_list_2.txt> ...\nInput lists are single-column gene IDs.")
}
# args <- c('out', 'cultivated', 'cordophanus', 'mucosospermus', 'amarus', 'colocynthis', 'three_wild');

out_pref <- as.character(args[1]);
args <- args[-1];

# Read gene lists and store in a named list
gene_lists <- lapply(args, read_gene_list)
names(gene_lists) <- args  # Keep file names as names

# Create a binary matrix for UpSet plot
all_genes <- unique(unlist(gene_lists))
binary_matrix <- as.data.frame(matrix(0, nrow = length(all_genes), ncol = length(gene_lists)))
colnames(binary_matrix) <- names(gene_lists)
rownames(binary_matrix) <- all_genes

for (i in seq_along(gene_lists)) {
  binary_matrix[rownames(binary_matrix) %in% gene_lists[[i]], i] <- 1
}

# Convert to data frame for UpSet
binary_matrix <- binary_matrix %>% rownames_to_column("Gene")

# Order columns by size
# set_size <- colSums(binary_matrix[, -1])
# binary_matrix <- binary_matrix[, c("Gene", names(sort(set_size, decreasing = TRUE))) ]

# Save plot to file
svg(paste0(out_pref, ".svg"), height=3, width= 6)
# upset(binary_matrix, sets = rev(colnames(binary_matrix)[-1]), keep.order=TRUE, order.by = "freq", nintersects= NA) # main.bar.color = "blue", sets.bar.color = "blue";
upset(binary_matrix, sets = rev(colnames(binary_matrix)[-1]), keep.order=TRUE, order.by = "freq", nintersects= 30) # main.bar.color = "blue", sets.bar.color = "blue";
dev.off()

