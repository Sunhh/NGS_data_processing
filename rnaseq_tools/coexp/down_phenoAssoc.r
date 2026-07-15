#!/usr/bin/env Rscript

# Associate WGCNA modules with sample traits, from a network built by run_wgcna.r.
# Merges the old down_phenoAssoc_pearson.r / down_phenoAssoc_binary.r (they differed only
# in --corType / --binaryTrait). Writes a labelled module-trait heatmap (pdf) and a
# module x trait correlation + p-value table (tsv).
#
# Build the network once (run_wgcna.r) then run this for each trait set — no rebuild.

suppressMessages({ library(optparse) })
opt_list <- list(
  make_option(c("-n","--network"), type="character", help="*_network_ForPhenoAssoc.RData from run_wgcna.r"),
  make_option(c("-p","--pheno"),   type="character", help="Trait table: row1 header, col1 = sample ID, other cols = traits (NA allowed)"),
  make_option(c("-o","--out"),     type="character", help="Output prefix (writes <out>.module_trait.heatmap.pdf and <out>.module_trait.tsv)"),
  make_option("--corType",     type="character", default="bicor", help="bicor | pearson [default: bicor]"),
  make_option("--binaryTrait", action="store_true", default=FALSE, help="Traits are binary (0/1); with bicor this sets robustY=FALSE"),
  make_option("--traitPval",   type="double", default=0.01, help="P-value threshold flagged 'Sig' in the heatmap [default: 0.01]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
for (r in c("network","pheno","out")) if (is.null(opt[[r]])) stop("Missing required --", r)
if (!opt$corType %in% c("bicor","pearson")) stop("--corType must be bicor|pearson")

suppressMessages({ library(WGCNA) }); options(stringsAsFactors = FALSE); enableWGCNAThreads()
robustY <- !(opt$corType == "bicor" && opt$binaryTrait)   # bicor + binary trait -> robustY FALSE

load(opt$network)   # MEs, SampleName, nSamples, moduleColors, ...
if (!exists("MEs")) stop("--network RData has no 'MEs' (was it made by run_wgcna.r?)")

allTraits <- read.table(opt$pheno, header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = "\t")
traitRows <- match(rownames(MEs), rownames(allTraits))
if (sum(!is.na(traitRows)) <= 1) stop("Too few samples shared between --network and --pheno")
datTraits <- allTraits[traitRows, , drop = FALSE]
rownames(datTraits) <- rownames(MEs)
good <- vapply(datTraits, function(x) sd(x, na.rm = TRUE) > 0, logical(1))
if (!any(good)) stop("No trait column has non-zero variance in the shared samples")
datTraits <- datTraits[, good, drop = FALSE]
message("[Msg] ", ncol(MEs), " modules x ", ncol(datTraits), " traits over ", nrow(MEs), " samples (corType=", opt$corType, ", binaryTrait=", opt$binaryTrait, ")")

if (opt$corType == "bicor") {
  s <- WGCNA::bicorAndPvalue(MEs, use = "pairwise.complete.obs")
  t <- WGCNA::bicorAndPvalue(MEs, datTraits, use = "pairwise.complete.obs", robustY = robustY)
  cor.self <- s$bicor; pval.self <- s$p; cor.trait <- t$bicor; pval.trait <- t$p
} else {
  cor.self  <- WGCNA::cor(MEs, use = "pairwise.complete.obs");            pval.self  <- WGCNA::corPvalueStudent(cor.self, nrow(MEs))
  cor.trait <- WGCNA::cor(MEs, datTraits, use = "pairwise.complete.obs"); pval.trait <- WGCNA::corPvalueStudent(cor.trait, nrow(MEs))
}

## ---- table: module x trait correlation + pvalue ------------------------------
tab <- data.frame(Module = rownames(cor.trait), check.names = FALSE)
for (tr in colnames(cor.trait)) { tab[[paste0("cor.", tr)]] <- cor.trait[, tr]; tab[[paste0("pval.", tr)]] <- pval.trait[, tr] }
write.table(tab, paste0(opt$out, ".module_trait.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- heatmap -----------------------------------------------------------------
mkText <- function(cor, pval) { tag <- ifelse(pval < opt$traitPval, "Sig", "")
  m <- paste0(signif(cor, 2), "\n(", signif(pval, 1), ")\n", tag); dim(m) <- dim(cor); m }
txt.trait <- mkText(cor.trait, pval.trait)
txt.self  <- mkText(cor.self,  pval.self); diag(txt.self) <- ""
pdf(paste0(opt$out, ".module_trait.heatmap.pdf"), width = 10, height = 15)
par(mar = c(6, 8.8, 3, 2.2))
WGCNA::labeledHeatmap(cor.trait, xLabels = colnames(datTraits), yLabels = names(MEs), ySymbols = names(MEs),
  colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = txt.trait, setStdMargins = FALSE, cex.text = 0.8, zlim = c(-1,1), main = "Module-trait relationships")
WGCNA::labeledHeatmap(cor.self, xLabels = names(MEs), yLabels = names(MEs), ySymbols = names(MEs),
  colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = txt.self, setStdMargins = FALSE, cex.text = 0.8, zlim = c(-1,1), main = "Module-module relationships")
dev.off()
message("[Rec] wrote ", opt$out, ".module_trait.{heatmap.pdf,tsv}")
