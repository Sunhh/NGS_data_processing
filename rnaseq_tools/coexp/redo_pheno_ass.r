#!/home/Sunhh/bin/Rscript
argvs <- commandArgs( trailingOnly=TRUE );
fn_RData <- as.character( argvs[1] );
fn_pheno <- as.character( argvs[2] );
fn_ModTrait_HmapPdf <- as.character( argvs[3] );
# fn_pheno <- '../input/dat2_pheno.add1'
# fn_ModTrait_HmapPdf <- 'dat2_module_trait_add1.heatmap.pdf'
### 'NA' in fn_pheno is accepted with WGCNA::cor(use="pairwise.complete.obs"), but this tends to provide lower p-values than using zero.

moduleTraitPval_thres <- 0.01

library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

load(fn_RData)
allTraits <- read.table( fn_pheno, header=T, row.names=1, stringsAsFactors=F, sep="\t" )
traitRows <- match(SampleName, rownames(allTraits));
datTraits <- allTraits[traitRows, ];
if ( ! is.data.frame(datTraits) ) {
  datTraits <- data.frame(datTraits)
  colnames(datTraits) <- colnames(allTraits)
}
rownames(datTraits) <- rownames(allTraits)[traitRows]
cor.mod_self              <-  WGCNA::cor(MEs, use = "pairwise.complete.obs")
pval.mod_self             <-  WGCNA::corPvalueStudent(cor.mod_self, nSamples)
sigTag.mod_self           <-  pval.mod_self < moduleTraitPval_thres
sigTag.mod_self[ sigTag.mod_self == FALSE ] <- ""
sigTag.mod_self[ sigTag.mod_self == TRUE  ] <- "Sig"
textMatrix.mod_self       <- paste(signif(cor.mod_self, 2), "\n(", signif(pval.mod_self, 1), ")\n", sigTag.mod_self, sep = "");
dim(textMatrix.mod_self)  <- dim(cor.mod_self)

cor.mod_trait             <- WGCNA::cor(MEs, datTraits, use = "pairwise.complete.obs")
pval.mod_trait            <- WGCNA::corPvalueStudent(cor.mod_trait, nSamples)
sigTag.mod_trait          <- pval.mod_trait < moduleTraitPval_thres
sigTag.mod_trait[ sigTag.mod_trait == FALSE ] <- ""
sigTag.mod_trait[ sigTag.mod_trait == TRUE  ] <- "Sig"
textMatrix.mod_trait      <- paste(signif(cor.mod_trait, 2), "\n(", signif(pval.mod_trait, 1), ")\n", sigTag.mod_self, sep = "")
dim(textMatrix.mod_trait) <- dim(cor.mod_trait)

pdf( fn_ModTrait_HmapPdf, width=10, height=15)
par(mar = c(6, 8.8, 3, 2.2))
WGCNA::labeledHeatmap(
  Matrix        = cor.mod_trait,             #相关系数
  xLabels       = colnames(datTraits),       #x轴为表型
  yLabels       = names(MEs),                #y轴为模块
  ySymbols      = names(MEs),
  colorLabels   = FALSE,
  colors        = blueWhiteRed(50),
  textMatrix    = textMatrix.mod_trait,      #每个单元格的内容
  setStdMargins = FALSE,
  cex.text      = 0.8,
  zlim          = c(-1,1),
  main          = "Module-trait relationships"
)
WGCNA::labeledHeatmap(
  Matrix        = cor.mod_self,             #相关系数
  xLabels       = names(MEs),               #x轴为表型
  yLabels       = names(MEs),               #y轴为模块
  ySymbols      = names(MEs),
  colorLabels   = FALSE,
  colors        = blueWhiteRed(50),
  textMatrix    = textMatrix.mod_self,      #每个单元格的内容
  setStdMargins = FALSE,
  cex.text      = 0.8,
  zlim          = c(-1,1),
  main          = "Module-module relationships"
)
dev.off()

