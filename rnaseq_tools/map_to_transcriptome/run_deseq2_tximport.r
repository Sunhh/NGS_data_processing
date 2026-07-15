#!/usr/bin/env Rscript

# tximport-based DESeq2, FDR only (one column per comparison), for the transcriptome
# pipeline. Uses DESeqDataSetFromTximport so the transcript-length information from
# Salmon/kallisto is respected. The up/down/not DEG decision is NOT made here; it is
# made downstream by DEG_byList_vClaude.pl (TPM fold change + adjustable thresholds).
#
# Output: GeneID + one padj column per comparison named ds.<group1>_VS_<group2>,
#         exactly the --fdr table shape DEG_byList_vClaude.pl consumes.

suppressMessages({ library(optparse); library(tximport); library(DESeq2) })

option_list <- list(
  make_option(c("-r","--rds"),         type="character", help="tximport RDS from merge_tx_quant_to_gene.r"),
  make_option(c("-m","--metadata"),    type="character", help="Sample table: columns 'sample' and 'group' (tab, header)"),
  make_option(c("-c","--compareList"), type="character", help="method<TAB>group1<TAB>group2<TAB>name (method col ignored; DESeq2 always). Same file as vClaude --pair"),
  make_option(c("-o","--output"),      type="character", help="Output FDR table")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (r in c("rds","metadata","compareList","output")) if (is.null(opt[[r]])) stop("Missing required --", r)

txi  <- readRDS(opt$rds)
meta <- read.table(opt$metadata, sep="\t", header=TRUE, stringsAsFactors=FALSE)
if (!all(c("sample","group") %in% colnames(meta))) stop("metadata needs 'sample' and 'group' columns")
rownames(meta) <- meta$sample
cmp  <- read.table(opt$compareList, sep="\t", header=FALSE, stringsAsFactors=FALSE, comment.char="")

allGenes <- rownames(txi$counts)
out <- data.frame(GeneID = allGenes, stringsAsFactors = FALSE, check.names = FALSE)

for (i in seq_len(nrow(cmp))) {
  g1 <- trimws(cmp[i,2]); g2 <- trimws(cmp[i,3])          # baseline, treatment
  colname <- paste0("ds.", g1, "_VS_", g2)
  use <- rownames(meta)[meta$group %in% c(g1, g2)]
  use <- intersect(colnames(txi$counts), use)
  sub <- meta[use, , drop=FALSE]
  n1 <- sum(sub$group == g1); n2 <- sum(sub$group == g2)
  if (n1 < 2 || n2 < 2) { message("[Wrn] <2 replicates for ", g1, "(", n1, ")/", g2, "(", n2, "); skip ", colname); next }
  txs <- txi
  txs$counts    <- txi$counts[,    use, drop=FALSE]
  txs$abundance <- txi$abundance[, use, drop=FALSE]
  txs$length    <- txi$length[,    use, drop=FALSE]
  sub$group <- factor(sub$group, levels = c(g1, g2))
  dds <- DESeqDataSetFromTximport(txs, colData = sub, design = ~ group)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", g2, g1))
  padj <- res$padj; names(padj) <- rownames(res)
  out[[colname]] <- padj[allGenes]
  message("[Msg] Done ", colname, " (", g1, " x", n1, " vs ", g2, " x", n2, ")")
}
if (ncol(out) < 2) stop("No comparison produced FDRs (all skipped?).")
write.table(out, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
