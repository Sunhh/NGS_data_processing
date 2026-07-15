#!/usr/bin/env Rscript

# Build a WGCNA co-expression network from an expression matrix and save the module
# assignment, module eigengenes (for trait association), per-gene module membership (KME),
# and diagnostic plots. Trait association is a SEPARATE step (down_phenoAssoc.r), so a
# network is built once and can be associated with several trait sets without rebuilding.
#
# Merges the old run_wgcna_signed.r / run_wgcna_unsigned.r (they differed only in
# --networkType) and applies current WGCNA best practice by default: build on log2
# expression (--transform) with the robust bicor correlation (--corType) on a signed
# network (--networkType). Use --transform none / --corType pearson / --networkType
# unsigned to reproduce older behaviour.

suppressMessages({ library(optparse) })

opt_list <- list(
  make_option(c("-e","--expr"),   type="character", help="Expression table: row1 header, col1 = gene ID, other cols = samples (NA as '-')"),
  make_option(c("-o","--outdir"), type="character", help="Output directory"),
  make_option(c("-p","--prefix"), type="character", help="Output file prefix"),
  make_option("--networkType",  type="character", default="signed",  help="signed | unsigned | 'signed hybrid' [default: signed]"),
  make_option("--corType",      type="character", default="bicor",   help="bicor | pearson [default: bicor]"),
  make_option("--transform",    type="character", default="log2",    help="Input transform before building the network: log2 | none [default: log2]"),
  make_option("--softPower",    type="integer",   default=NA,        help="Soft-thresholding power; omitted = pick automatically"),
  make_option("--minModuleSize",type="integer",   default=30L,       help="Minimum genes per module [default: 30]"),
  make_option("--mergeCutHeight",type="double",   default=0.25,      help="Module-merge cut height; larger = fewer modules [default: 0.25]"),
  make_option("--maxBlockSize", type="integer",   default=50000L,    help="Max block size; keep >= gene number to use one block [default: 50000]"),
  make_option("--threads",      type="integer",   default=0L,        help="WGCNA threads; 0 = auto [default: 0]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
for (r in c("expr","outdir","prefix")) if (is.null(opt[[r]])) stop("Missing required --", r)
if (!opt$networkType %in% c("signed","unsigned","signed hybrid")) stop("--networkType must be signed|unsigned|'signed hybrid'")
if (!opt$corType %in% c("bicor","pearson")) stop("--corType must be bicor|pearson")
if (!opt$transform %in% c("log2","none")) stop("--transform must be log2|none")

suppressMessages({ library(WGCNA) })
options(stringsAsFactors = FALSE)
if (opt$threads > 0) enableWGCNAThreads(opt$threads) else enableWGCNAThreads()

corType      <- opt$corType
maxPOutliers <- ifelse(corType == "bicor", 0.05, 1)
corFnc       <- ifelse(corType == "bicor", "bicor", "cor")
networkType  <- opt$networkType
RsquaredCut  <- c(0.9, 0.8)
powers       <- if (networkType %in% c("unsigned","signed hybrid")) 1:15 else 1:30

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
oo <- file.path(opt$outdir, opt$prefix)

## ---- load, transform, filter -------------------------------------------------
expr0 <- read.table(opt$expr, header = TRUE, row.names = 1, sep = "\t", na.strings = "-")
expr0 <- if (opt$transform == "log2") log2(expr0 + 0.01) else expr0
datExprOri <- t(expr0)
nGenes0 <- ncol(datExprOri)
keepSample <- apply(is.na(data.frame(datExprOri)), 1, sum) < 0.1 * nGenes0
datExpr    <- datExprOri[keepSample, ]
SampleName <- rownames(datExprOri)[keepSample]
nSamples   <- nrow(datExpr)

vExpr  <- apply(as.matrix(datExpr), 2, var, na.rm = TRUE)
nMiss  <- apply(is.na(as.matrix(datExpr)), 2, sum)
goodEx <- apply(as.matrix(datExpr), 2, max, na.rm = TRUE) >= ifelse(opt$transform == "log2", log2(1 + 0.01), 1)
keepGenes <- vExpr > 0 & nMiss < 0.1 * nSamples & goodEx
datExpr   <- datExpr[, keepGenes]
message("[Msg] kept ", nSamples, " samples x ", ncol(datExpr), " genes (from ", nGenes0, ")")
save(datExpr, file = paste0(oo, "_datExprInput.RData"))

## ---- sample clustering -------------------------------------------------------
pdf(paste0(oo, "_sample_hclust.pdf"), width = 7, height = 5)
plot(hclust(dist(datExpr), method = "average"), xlab = "", sub = "", cex = 0.7, main = "Sample clustering", hang = -1)
dev.off()

## ---- soft-threshold power ----------------------------------------------------
sft_power <- opt$softPower
if (is.na(sft_power)) {
  sft <- pickSoftThreshold(datExpr, RsquaredCut = RsquaredCut[1], powerVector = powers, networkType = networkType, verbose = 0)
  if (!is.na(sft$powerEstimate)) {
    sft_power <- sft$powerEstimate
  } else {
    for (rc in RsquaredCut[-1]) {
      hit <- which(sft$fitIndices$SFT.R.sq >= rc)
      if (length(hit) > 0) { sft_power <- sft$fitIndices$Power[hit[1]]; message("[Msg] soft power ", sft_power, " by R^2 >= ", rc); break }
    }
    if (is.na(sft_power)) {   # empirical fallback by sample size
      base <- if (networkType %in% c("unsigned","signed hybrid")) c(9,8,7,6) else c(18,16,14,12)
      sft_power <- if (nSamples < 20) base[1] else if (nSamples < 30) base[2] else if (nSamples < 40) base[3] else base[4]
      message("[Msg] soft power ", sft_power, " by empirical value (nSamples=", nSamples, ", ", networkType, ")")
    }
  }
  pdf(paste0(oo, "_sftThreshold.pdf"), width = 7, height = 4)
  par(mfrow = c(1,2)); cex1 <- 0.9
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft power", ylab = "Scale-free R^2", type = "n", main = "Scale independence", ylim = c(-1,1))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = cex1, col = "red"); abline(h = 0.9, col = "red"); abline(h = 0.8, col = "blue")
  plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft power", ylab = "Mean connectivity", type = "n", main = "Mean connectivity")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")
  dev.off()
}
message("[Msg] soft power used: ", sft_power)

## ---- build network -----------------------------------------------------------
net <- blockwiseModules(
  datExpr, power = sft_power, maxBlockSize = opt$maxBlockSize,
  TOMType = "unsigned", networkType = networkType, deepSplit = 2,
  minModuleSize = opt$minModuleSize, mergeCutHeight = opt$mergeCutHeight,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE, saveTOMFileBase = paste0(oo, "_TOM"), loadTOMs = TRUE,
  corType = corType, maxPOutliers = maxPOutliers, verbose = 0
)

moduleColors <- labels2colors(net$colors)
moduleLabels <- net$colors
MEs          <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
rownames(MEs) <- SampleName
net.MEs      <- net$MEs
net.geneTree <- net$dendrograms[[1]]
save(moduleColors, moduleLabels, MEs, net.MEs, net.geneTree, SampleName, nSamples,
     file = paste0(oo, "_network_ForPhenoAssoc.RData"))
message("[Msg] ", length(unique(moduleColors)), " modules (incl. grey)")

## ---- KME (module membership) table ------------------------------------------
kme  <- signedKME(datExpr, MEs, corFnc = corFnc, corOptions = "use = 'p'", outputColumnName = "ME")
midx <- match(paste0("ME", moduleColors), colnames(kme))
kmeInM <- vapply(seq_len(nrow(kme)), function(i) kme[i, midx[i]], numeric(1))
colnames(kme) <- paste0("gMM_", substring(colnames(kme), 3))
conn <- intramodularConnectivity.fromExpr(datExpr, moduleColors, corFnc = corFnc, corOptions = "use = 'p'", networkType = networkType)
kme_tbl <- data.frame(EleID = colnames(datExpr), KME = kmeInM, ModuleID = moduleColors,
                      absKME = abs(kmeInM), ModuleNum = net$colors, conn[, c("kWithin","kOut")], kme)
kme_tbl <- kme_tbl[order(kme_tbl$ModuleNum, -kme_tbl$absKME), ]
write.table(kme_tbl, file = paste0(oo, "_KME.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## ---- diagnostic plots (pdf; peripheral ones are best-effort) -----------------
expColor <- t(numbers2colors(datExpr, colors = blueWhiteRed(100), naColor = "grey"))
colnames(expColor) <- rownames(datExpr)
pdf(paste0(oo, "_dendroColors.pdf"), height = 14, width = 14)
plotDendroAndColors(net$dendrograms[[1]], colors = cbind(moduleColors[net$blockGenes[[1]]], expColor),
                    c("Module", colnames(expColor)), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, cex.rowText = 0.5)
dev.off()
pdf(paste0(oo, "_adjHeatmap.pdf"), height = 7, width = 7)
plotEigengeneNetworks(MEs, "Eigengene adjacency", plotDendrograms = TRUE, marDendro = c(4,4,2,4), marHeatmap = c(4,10,2,10))
dev.off()

tryCatch({
  d <- paste0(oo, "_ME_barplot"); dir.create(d, showWarnings = FALSE)
  for (module in substring(colnames(MEs), 3)) {
    if (module == "grey") next
    pdf(file.path(d, paste0(module, ".express.barplot.pdf")), height = 7, width = 9)
    par(mfrow = c(2,1), mar = c(0.3, 6.5, 5, 3))
    plotMat(t(scale(datExpr[, moduleColors == module])), rlabels = FALSE, main = module, cex.main = 2, clabels = SampleName)
    par(mar = c(5, 4.2, 0, 0.7))
    barplot(MEs[, paste0("ME", module)], col = module, main = "", cex.main = 2, ylab = "eigengene expression", xlab = "sample")
    dev.off()
  }
}, error = function(e) message("[Wrn] ME barplots skipped: ", conditionMessage(e)))

tryCatch({
  d <- paste0(oo, "_cytoscape/"); dir.create(d, showWarnings = FALSE)
  load(paste0(oo, "_TOM-block.1.RData"))   # loads TOM
  TOM.mat <- as.matrix(TOM)
  for (module in substring(colnames(MEs), 3)) {
    if (module == "grey") next
    inM <- moduleColors == module
    modProbes <- colnames(datExpr)[inM]
    modTOM <- TOM.mat[inM, inM]; dimnames(modTOM) <- list(modProbes, modProbes)
    exportNetworkToCytoscape(modTOM, edgeFile = paste0(d, "input-edges-", module, ".txt"),
      nodeFile = paste0(d, "input-nodes-", module, ".txt"), weighted = TRUE, threshold = 0.02,
      nodeNames = modProbes, nodeAttr = moduleColors[inM])
  }
}, error = function(e) message("[Wrn] cytoscape export skipped: ", conditionMessage(e)))

message("[Rec] network built. Associate modules with traits via down_phenoAssoc.r using ",
        oo, "_network_ForPhenoAssoc.RData")
