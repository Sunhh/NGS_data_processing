# Load required packages
suppressMessages({library(optparse)})

# Define command-line arguments
option_list <- list(
  make_option(c("-o", "--go_obo"), type="character", default=NULL, 
              help="Path to GO OBO file", metavar="character"),
  make_option(c("-p", "--output_prefix"), type="character", default="out", 
              help="Prefix for output files", metavar="character")
)

three_root_GOs <- c('GO:0008150', 'GO:0003674', 'GO:0005575'); # BP, MF, CC;

# Parse command-line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
# opt <- list(
#  go_obo = '20250108-go.obo',
#  output_prefix = 'to1'
# );


# Check if required arguments are provided
if (is.null(opt$go_obo)) {
  print_help(opt_parser)
  stop("All --go_obo arguments must be provided.", call.=FALSE);
}

if (is.null(opt$output_prefix)) {
  opt$output_prefix <- 'out'
}

suppressMessages({
library(ontologyIndex)
library(dplyr)
});

# Load GO OBO file
obo <- get_ontology(opt$go_obo, extract_tags = "everything", propagate_relationships= c('is_a', 'part_of'))

# Make GO mapping table.
#   There can be 3%-6% differences between GO enrichment sets from the OBO file-based R program and Blast2GO software.
#   The EBI QuickGO supports the set from OBO file-based R program.
### Map GO IDs with alternative GO IDs (combined), root, GO name, GO description, ancestors.
good_GO_IDs <- obo$id[grepl("^GO:", obo$id)]
obo_df1 <- do.call(rbind, lapply(good_GO_IDs, function(term_id) {
  # Get alternative GO IDs
  alt_ids <- obo$alt_id[[term_id]]
  if (is.null(alt_ids) || length(alt_ids) == 0) {
    alt_ids <- NA
  }
  # Get root ontology;
  rootGO <- obo$namespace[[term_id]];
  # Get GO name
  nameGO <- obo$name[[term_id]]
  # Get GO description. # Some of the "def" extractions are shorter than correct. It should be OK for large-scale analysis.
  descGO <- obo$def[[term_id]]
  # Get ancestors; (including self without roots)
  ancestors <- obo$ancestors[[term_id]];
  ancestors <- ancestors[ !(ancestors %in% three_root_GOs) ];
  if (length(ancestors) == 0) { ancestors <- c(term_id) }
  data.frame(
    GO_ID = term_id,
    ALT_ID = paste(alt_ids, collapse=";"),
    Root_GO = rootGO,
    GO_name  = nameGO,
    GO_desc  = descGO, 
    Ancestors= paste(ancestors, collapse=";"),
    stringsAsFactors = FALSE)
}))

# Rename ontology levels
obo_df1$Root_GO <- recode(obo_df1$Root_GO,
                                "biological_process" = "Biological Process",
                                "molecular_function" = "Molecular Function",
                                "cellular_component" = "Cellular Component")

# Output tables.
### Merged table for all information.
write.table(obo_df1, file= paste0(opt$output_prefix, '.tab'), sep= "\t", quote= FALSE, row.names= FALSE, col.names=TRUE);

