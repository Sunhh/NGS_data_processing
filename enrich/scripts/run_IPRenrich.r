# Load required packages
suppressMessages({library(optparse)})

# Define command-line arguments
option_list <- list(
  make_option(c("-g", "--gene_list"), type="character", default=NULL, 
              help="Path to gene list file", metavar="character"),
  make_option(c("-a", "--ipr_annot_3col"), type="character", default=NULL, 
              help="Path to GO annotation file with 4 columns: geneID, GO ID, root GO, GO name", metavar="character"),
  make_option(c("-p", "--output_prefix"), type="character", default="out", 
              help="Prefix for output files", metavar="character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


# Check if required arguments are provided
if (is.null(opt$gene_list) || is.null(opt$ipr_annot_3col)) {
  print_help(opt_parser)
  stop("All --gene_list, --ipr_annot_3col arguments must be provided.", call.=FALSE)
}

if (is.null(opt$output_prefix)) {
  opt$output_prefix <- 'out'
}



suppressMessages({
library(clusterProfiler)
library(ggplot2)
library(GO.db)
library(ontologyIndex)
library(dplyr)
});

# Read input files
gene_list <- read.table(opt$gene_list, header=FALSE, stringsAsFactors=FALSE)[,1]  # Your gene IDs
go_annotations <- read.table(opt$ipr_annot_3col, header=FALSE, stringsAsFactors=FALSE, quote= "", fill= TRUE, sep="\t", blank.lines.skip = FALSE); # GO annotations
colnames(go_annotations)[1:3] <- c('Gene', 'IPR_ID', 'IPR_name');

# Generate input variables for enricher().
in_uniqGOs <- go_annotations %>%
 select(-Gene) %>%
 distinct(IPR_ID, .keep_all = TRUE)
### TERM2NAME: in_uniqGOs     ==>> IPR_ID to IPR_name;
### TERM2GENE: go_annotations ==>> IPR_ID to Gene;    # which includes non-annotated gene IDs.

# Perform GO enrichment analysis using custom GO annotation file
ego <- enricher(gene_list, TERM2GENE = go_annotations[, c('IPR_ID', 'Gene')], TERM2NAME = in_uniqGOs[,c('IPR_ID', 'IPR_name')], pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize= 3)

# Merge results with GO categories
ego@result <- merge(ego@result, in_uniqGOs, by.x="ID", by.y="IPR_ID", all.x=TRUE, sort= FALSE) %>% select(-IPR_name)

# Calculate new GeneRatio as (GO genes in input) / (GO genes in background)
ego@result <- ego@result %>%
  mutate(GeneRatio_in_term = paste0(Count, "/", as.numeric(sub("/.*", "", BgRatio))))

# Save the results table
output_table <- paste0(opt$output_prefix, ".tsv")
write.table(ego@result, file=output_table, sep="\t", quote=FALSE, row.names=FALSE)

if (sum(ego@result$p.adjust < 0.05) > 0) {
# Filter results for p.adjust < 0.05 and sort p.adjust
filtered_results <- ego@result %>% filter(p.adjust < 0.05) %>%
  mutate(GeneRatio = GeneRatio_in_term)

plot_ego <- ego
plot_ego@result <- filtered_results

# Save the plot to an SVG file
output_svg <- paste0(opt$output_prefix, ".svg")
p1 <- dotplot(plot_ego, showCategory=30) +
 theme(panel.grid = element_blank())
svg(output_svg, width=8, height=5)
print(p1)
dev.off()

if (sum(filtered_results$Count >= 10) > 0) {
filtered_res2 <- filtered_results %>% filter(Count >= 10)
plot_ego@result <- filtered_res2
p1 <- dotplot(plot_ego, showCategory=30) +
 theme(panel.grid = element_blank())
output_svg2 <- paste0(opt$output_prefix, "-testLg10.svg");
svg(output_svg2, width=8, height=5)
print(p1)
dev.off()
}
}

