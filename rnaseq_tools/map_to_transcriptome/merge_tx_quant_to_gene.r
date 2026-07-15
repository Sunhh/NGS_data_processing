#!/usr/bin/env Rscript

# Merge per-sample transcript quantifications to gene level with tximport.
# Works for both Salmon and kallisto (choose with --type); the rest of the
# transcriptome DEG chain (run_deseq2_tximport.r -> DEG_byList_vClaude.pl) is
# identical whichever quantifier produced the per-sample directories.
#
# Outputs (prefix from --out):
#   <out>.gene_counts.tsv    gene x sample, lengthScaledTPM counts (for DESeq2)
#   <out>.gene_tpm.tsv       gene x sample, TPM (fed to vClaude as --tpm)
#   <out>.gene_lengths.tsv   gene x sample, average transcript length
#   <out>.RDS                the tximport object (input to run_deseq2_tximport.r)

suppressMessages({ library(optparse); library(tximport); library(readr) })

option_list <- list(
  make_option(c("-s","--samples"),  type="character", help="File listing per-sample quant directories, one per line"),
  make_option(c("-g","--tx2gene"),  type="character", help="transcript_id<TAB>gene_id table (with header row)"),
  make_option(c("-o","--out"),      type="character", help="Output prefix"),
  make_option(c("-t","--type"),     type="character", default="salmon", help="Quantifier: salmon | kallisto [default: salmon]"),
  make_option(c("--quant_file"),    type="character", default=NULL, help="Quant filename inside each dir [salmon: quant.sf; kallisto: abundance.h5]"),
  make_option(c("--countsFromAbundance"), type="character", default="lengthScaledTPM",
              help="tximport countsFromAbundance: no|scaledTPM|lengthScaledTPM [default: lengthScaledTPM]")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (r in c("samples","tx2gene","out")) if (is.null(opt[[r]])) stop("Missing required --", r)
opt$type <- tolower(opt$type)
if (!opt$type %in% c("salmon","kallisto")) stop("--type must be 'salmon' or 'kallisto'")

quant_file <- opt$quant_file
if (is.null(quant_file)) quant_file <- if (opt$type == "salmon") "quant.sf" else "abundance.h5"

sample_dirs <- readLines(opt$samples)
sample_dirs <- sample_dirs[nchar(trimws(sample_dirs)) > 0]
files <- file.path(sample_dirs, quant_file)
names(files) <- basename(sample_dirs)
missing <- !file.exists(files)
if (any(missing)) stop("Missing ", quant_file, " in:\n", paste(sample_dirs[missing], collapse="\n"))

tx2gene <- read_tsv(opt$tx2gene, col_names = TRUE, show_col_types = FALSE)

txi <- tximport(files, type = opt$type, tx2gene = tx2gene,
                countsFromAbundance = opt$countsFromAbundance)

wtab <- function(m, suf) write.table(data.frame(GeneID = rownames(m), m, check.names = FALSE),
                                     file = paste0(opt$out, suf), sep = "\t", quote = FALSE, row.names = FALSE)
wtab(txi$counts,    ".gene_counts.tsv")
wtab(txi$abundance, ".gene_tpm.tsv")
wtab(txi$length,    ".gene_lengths.tsv")
saveRDS(txi, paste0(opt$out, ".RDS"))
message("[Msg] tximport(", opt$type, ") merged ", ncol(txi$counts), " samples, ",
        nrow(txi$counts), " genes -> ", opt$out, ".{gene_counts,gene_tpm,gene_lengths}.tsv + .RDS")
