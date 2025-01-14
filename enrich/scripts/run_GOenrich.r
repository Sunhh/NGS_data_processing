# Load required packages
suppressMessages({library(optparse)})

# Define command-line arguments
option_list <- list(
  make_option(c("-g", "--gene_list"), type="character", default=NULL, 
              help="Path to gene list file", metavar="character"),
  make_option(c("-a", "--go_annot_4col"), type="character", default=NULL, 
              help="Path to GO annotation file with 4 columns: geneID, GO ID, root GO, GO name", metavar="character"),
  make_option(c("-p", "--output_prefix"), type="character", default="out", 
              help="Prefix for output files", metavar="character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# opt <- list(
#  gene_list = 'A_CLV-P_CLC.txt',
#  go_annot_4col = 'synFam.b2g.annot-GOinEnrich',
#  out_prefix = 'to1'
# );

# Check if required arguments are provided
if (is.null(opt$gene_list) || is.null(opt$go_annot_4col)) {
  print_help(opt_parser)
  stop("All --gene_list, --go_annot_4col arguments must be provided.", call.=FALSE)
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
go_annotations <- read.table(opt$go_annot_4col, header=FALSE, stringsAsFactors=FALSE, quote= "", fill= TRUE, sep="\t", blank.lines.skip = FALSE); # GO annotations
colnames(go_annotations)[1:4] <- c('Gene', 'GO_ID', 'GO_root', 'GO_name');

# Generate input variables for enricher().
in_uniqGOs <- go_annotations %>%
 select(-Gene) %>%
 distinct(GO_ID, .keep_all = TRUE)
### TERM2NAME: in_uniqGOs     ==>> GO_ID to GO_name;
### TERM2GENE: go_annotations ==>> GO_ID to Gene;    # which includes non-annotated gene IDs.

# Perform GO enrichment analysis using custom GO annotation file
ego <- enricher(gene_list, TERM2GENE = go_annotations[, c('GO_ID', 'Gene')], TERM2NAME = in_uniqGOs[,c('GO_ID', 'GO_name')], pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, minGSSize= 3)

# Merge results with GO categories
ego@result <- merge(ego@result, in_uniqGOs, by.x="ID", by.y="GO_ID", all.x=TRUE, sort= FALSE)

# Calculate new GeneRatio as (GO genes in input) / (GO genes in background)
ego@result <- ego@result %>%
  mutate(GeneRatio_in_term = paste0(Count, "/", as.numeric(sub("/.*", "", BgRatio))))

# Save the results table
output_table <- paste0(opt$output_prefix, ".tsv")
output_Tag <- ifelse(ego@result$p.adjust < 0.05, "OVER", "")
write.table(data.frame(Tag=output_Tag, ego@result), file=output_table, sep="\t", quote=FALSE, row.names=FALSE)

if (sum(ego@result$p.adjust < 0.05) > 0) {
# Filter results for p.adjust < 0.05 and sort by Ontology and p.adjust
filtered_results <- ego@result %>% filter(p.adjust < 0.05) %>%
  mutate(GeneRatio = GeneRatio_in_term, Ontology = GO_root)

# Plot dot plot with facets
filtered_results$Ontology <- factor(filtered_results$Ontology, levels = c("Biological Process", "Molecular Function", "Cellular Component"))

plot_ego <- ego
plot_ego@result <- filtered_results

# Save the plot to an SVG file
p1 <- dotplot(plot_ego, showCategory=30, split = 'Ontology') +
 facet_grid(Ontology~., scale="free", space= 'free_y') +
 theme(panel.grid = element_blank())
output_svg <- paste0(opt$output_prefix, ".svg")
svg(output_svg, width=10, height=10)
print(p1);
dev.off()

if (sum(filtered_results$Count >= 10) > 0) {
filtered_res2 <- filtered_results %>% filter(Count >= 10)
plot_ego@result <- filtered_res2
p1 <- dotplot(plot_ego, showCategory=30, split = 'Ontology') +
 facet_grid(Ontology~., scale="free", space= 'free_y') +
 theme(panel.grid = element_blank())
output_svg2 <- paste0(opt$output_prefix, "-testLg10.svg");
svg(output_svg2, width=10, height=10)
print(p1)
dev.off()
}
}

