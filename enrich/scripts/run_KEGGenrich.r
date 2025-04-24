#!/usr/bin/env Rscript
suppressMessages({library(optparse)})
# Define command-line arguments
option_list <- list(
  make_option(c("-g", "--gene_list"), type="character", default=NULL,
              help="Path to gene list file", metavar="character"),
  make_option(c("-a", "--kegg_annot_2col"), type="character", default=NULL,
              help="Path to KEGG annotation file with 2 columns: geneID, KEGG ID", metavar="character"),
  make_option(c("-p", "--output_prefix"), type="character", default="out",
              help="Prefix for output files", metavar="character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$gene_list) || is.null(opt$kegg_annot_2col)) {
  print_help(opt_parser)
  stop("All --gene_list, --kegg_annot_2col arguments must be provided.", call.=FALSE)
}

if (is.null(opt$output_prefix)) {
  opt$output_prefix <- 'out'
}

suppressMessages({
library(KEGGREST)
library(clusterProfiler)
library(dplyr)
})


# Load KEGG annotation (KO).
#   k_tbl   <- read.table("97103_v3-ghostKID.tsv", header = FALSE, col.names = c("gene","ko"), sep = "\t", stringsAsFactors = FALSE)
#   test_gene <- scan("RNAi_DEG", what= "")
test_gene <- scan(opt$gene_list, what= "")
k_tbl   <- read.table(opt$kegg_annot_2col, header = FALSE, col.names = c("gene","ko"), sep = "\t", stringsAsFactors = FALSE)

all_gene <- unique(k_tbl$gene)
k_tbl_df   <- subset(k_tbl, ko != "")


# Map KO to KEGG pathway.
ko2path <- keggLink("pathway", "ko") # Get information from KEGG server.
ko_only <- ko2path[ grepl("^path:ko", ko2path)   ] # Only keep KO pathways.
map_only <- ko2path[ grepl("^path:map", ko2path) ] # Only keep reference MAPs.
ko_ids <- sub("path:", "", unique(ko_only)) # Get all pathway IDs.
ko_desc <- keggList("pathway", "ko") # Get all pathway descriptions.
ko_desc_df <- data.frame(
  ko   = names(ko_desc),
  desc = ko_desc
)
path2k_df <- data.frame(
  pathway = sub("path:","", ko_only),
  ko      = sub("ko:","", names(ko_only))
)


# Generate gene to pathway annotation.
term2gene <- merge(k_tbl_df, path2k_df, by="ko")[,c("pathway","gene")]
term2name <- data.frame(ko=ko_ids) %>%
  left_join(ko_desc_df, by= "ko") %>%
  select(ko, desc) %>%
  arrange(ko, desc)
### gene2term <- k_tbl %>% select(gene) %>%
###   left_join(term2gene, by= "gene") %>%
###   select(gene, pathway) %>%
###   arrange(gene, pathway)

# Perform enrichment.
ego <- enricher(test_gene, TERM2GENE = term2gene, TERM2NAME = term2name, universe=all_gene, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize = 3)

# Calculate new GeneRatio as (GO genes in input) / (GO genes in background)
ego@result <- ego@result %>%
  mutate(GeneRatio_in_term = paste0(Count, "/", as.numeric(sub("/.*", "", BgRatio))))

# Save the results table
output_table <- paste0(opt$output_prefix, ".tsv")
output_Tag <- ifelse(ego@result$p.adjust < 0.05, "OVER", "")
write.table(data.frame(Tag=output_Tag, ego@result), file=output_table, sep="\t", quote=FALSE, row.names=FALSE)

