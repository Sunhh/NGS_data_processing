#!/usr/bin/env Rscript
suppressMessages({
library(VennDiagram)
library(tidyverse)
library(grid)
library(gridGraphics)
})

# Function to read gene lists from files
gene_list_from_file <- function(filename) {
  read_lines(filename) %>% unique()
}

base_colors <- c("blue", "red", "green", "yellow", "purple");

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("\nUsage: Rscript venn_plot.R out_prefix file1.txt file2.txt ... fileN.txt\n -- Got out_prefix.pdf plot.\n\n")
}

# Extract output prefix and file list
out_prefix <- args[1]
gene_files <- args[-1]

# Read gene lists
gene_lists <- lapply(gene_files, gene_list_from_file)
names(gene_lists) <- basename(gene_files)

# Generate Venn diagram and save as SVG
venn.plot <- venn.diagram(
  x = gene_lists,
  category.names = names(gene_lists),
  filename = NULL,
  disable.logging = TRUE,
  output = TRUE,
  fill = base_colors[1:length(gene_lists)],
  alpha = 0.5,
  cat.cex = 2,
  cex = 2
)

pdf(file = paste0(out_prefix, ".pdf"))
grid.newpage()
grid.draw(venn.plot)
dev.off()

#svg(filename = paste0(out_prefix, ".svg"))
#grid.newpage()
#grid.draw(venn.plot)
#dev.off()
