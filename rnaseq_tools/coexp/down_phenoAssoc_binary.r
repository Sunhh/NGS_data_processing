#!/home/Sunhh/bin/Rscript
# Values to be set before each run. 
# http://www.bio-info-trainee.com/2535.html
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html

argvs <- commandArgs( trailingOnly=TRUE );
if ( is.na(argvs[1]) ) {
  message("Rscript this.R   <in_network_ForPhenoAssoc.RData>   <input_pheno>   <out_module_trait.heatmap.pdf>   [moduleTraitPval_thres]");
  q(); 
}
fn_networkR                 <- as.character(argvs[1])
fn_pheno                    <- as.character(argvs[2])
fn_ModTrait_HmapPdf         <- as.character(argvs[3])
if ( is.na(argvs[4]) ) {
	moduleTraitPval_thres <- 0.01 ; 
} else {
	moduleTraitPval_thres <- as.character(argvs[4])
}

v_corType      <- 'bicor' ; # Could be 'pearson' or 'bicor';
v_binaryTrait  <- TRUE ; 
v_maxPOutliers <- ifelse( v_corType == 'bicor', 0.05, 1 ); 
v_robustY      <- ifelse( v_corType == 'bicor' && v_binaryTrait == TRUE, FALSE, TRUE ); 
# Problem parameters: corFnc ; maxPOutliers ; 
# 


### Start working; 
library(WGCNA)
enableWGCNAThreads() 

# Association : Here I want to add each sample as a pheno. 
### 'NA' in fn_pheno is accepted with WGCNA::cor(use="pairwise.complete.obs"), but this tends to provide lower p-values than using zero. 
### Varaibles required from previous analysis : 
##### v_corType : fixed. 
##### fn_pheno 
##### fn_ModTrait_HmapPdf
##### moduleTraitPval_thres
##### MEs     : Get from paste0(oopref, "_network_ForPhenoAssoc.RData"); 
##### SampleName # Should be rownames(MEs) or rownames(datExpr)
##### nSamples   # Should be nrow(MEs) or nrow(datExpr)
### So I want to save data for association anaylysis. 
load(fn_networkR); 
doAssoc <- function(
                    v_corType_Ass             = v_corType, 
		    fn_pheno_Ass              = fn_pheno, 
		    MEs_Ass                   = MEs, 
		    fn_ModTrait_HmapPdf_Ass   = fn_ModTrait_HmapPdf, 
		    moduleTraitPval_thres_Ass = moduleTraitPval_thres ) {
  allTraits <- read.table( fn_pheno_Ass, header=T, row.names=1, stringsAsFactors=F, sep="\t" )
  traitRows <- match(rownames(MEs_Ass), rownames(allTraits));
  # traitRows <- match(SampleName, rownames(allTraits));
  if (length(traitRows) <= 1) {
	  message("The available sample number is too low as ", length(traitRows), "\n"); 
	  return(1); 
  }
  datTraits.1 <- allTraits[traitRows, ];    
  if ( ! is.data.frame(datTraits.1) ) {
    datTraits.1 <- data.frame(datTraits.1)
    colnames(datTraits.1) <- colnames(allTraits)
  }
  rownames(datTraits.1) <- rownames(allTraits)[traitRows]
  datTraits.1.goodCols <- apply(datTraits.1, MARGIN=2, FUN=function(x){ sd(x, na.rm=T) > 0 } )
  datTraits  <-  datTraits.1[ , datTraits.1.goodCols]
  if ( ! is.data.frame(datTraits) ) {
    datTraits.1 <- data.frame(datTraits)
    colnames(datTraits) <- colnames(datTraits.1)[datTraits.1.goodCols]
  }

  if ( v_corType_Ass == 'bicor' ) {
    corP.mod_self              <- WGCNA::bicorAndPvalue(MEs_Ass, use = "pairwise.complete.obs")
    cor.mod_self               <- corP.mod_self$bicor
    pval.mod_self              <- corP.mod_self$p
    corP.mod_trait             <- WGCNA::bicorAndPvalue(MEs_Ass, datTraits, use = "pairwise.complete.obs", robustY = v_robustY)
    cor.mod_trait              <- corP.mod_trait$bicor
    pval.mod_trait             <- corP.mod_trait$p
  } else {
    cor.mod_self              <-  WGCNA::cor(MEs_Ass, use = "pairwise.complete.obs")
    pval.mod_self             <-  WGCNA::corPvalueStudent(cor.mod_self, nrow(MEs_Ass))
    cor.mod_trait             <-  WGCNA::cor(MEs_Ass, datTraits, use = "pairwise.complete.obs")
    pval.mod_trait            <-  WGCNA::corPvalueStudent(cor.mod_trait, nrow(MEs_Ass))
  }

  sigTag.mod_self           <-  pval.mod_self < moduleTraitPval_thres_Ass
  sigTag.mod_self[ sigTag.mod_self == FALSE ] <- ""
  sigTag.mod_self[ sigTag.mod_self == TRUE  ] <- "Sig"
  textMatrix.mod_self       <- paste(signif(cor.mod_self, 2), "\n(", signif(pval.mod_self, 1), ")\n", sigTag.mod_self, sep = "");
  dim(textMatrix.mod_self)  <- dim(cor.mod_self)
  for (i1 in 1:nrow(textMatrix.mod_self)) {
    textMatrix.mod_self[i1,i1] <- ''; 
  }

  sigTag.mod_trait          <- pval.mod_trait < moduleTraitPval_thres_Ass
  sigTag.mod_trait[ sigTag.mod_trait == FALSE ] <- ""
  sigTag.mod_trait[ sigTag.mod_trait == TRUE  ] <- "Sig"
  textMatrix.mod_trait      <- paste(signif(cor.mod_trait, 2), "\n(", signif(pval.mod_trait, 1), ")\n", sigTag.mod_trait, sep = "")
  dim(textMatrix.mod_trait) <- dim(cor.mod_trait)

  pdf( fn_ModTrait_HmapPdf_Ass, width=10, height=15)
  par(mar = c(6, 8.8, 3, 2.2))
  WGCNA::labeledHeatmap(
    Matrix        = cor.mod_trait,             #相关系数
    xLabels       = colnames(datTraits),       #x轴为表型
    yLabels       = names(MEs_Ass),                #y轴为模块
    ySymbols      = names(MEs_Ass),
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
    xLabels       = names(MEs_Ass),               #x轴为表型
    yLabels       = names(MEs_Ass),               #y轴为模块
    ySymbols      = names(MEs_Ass),
    colorLabels   = FALSE,
    colors        = blueWhiteRed(50),
    textMatrix    = textMatrix.mod_self,      #每个单元格的内容
    setStdMargins = FALSE,
    cex.text      = 0.8,
    zlim          = c(-1,1),
    main          = "Module-module relationships"
  )
  dev.off()

 #  return(); 
}# End of doAssoc() 

doAssoc(); 

message("All pheno associations done"); 

