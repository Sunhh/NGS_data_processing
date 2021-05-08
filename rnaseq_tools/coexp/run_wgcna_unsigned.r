#!/home/Sunhh/bin/Rscript
# Values to be set before each run. 
# http://www.bio-info-trainee.com/2535.html
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
fn_expr  <- 'input/dat1_rpkmMean'
fn_pheno <- 'input/dat1_pheno'
odir    <- 'wgcna_dat1'
opref   <- 'dat1'
sft_power <- NULL ;              # NULL means I'll estimate soft power threshold automatically. 

argvs <- commandArgs( trailingOnly=TRUE );
if ( is.na(argvs[1]) ) {
  message("Rscript this.R   <input_rpkmMean>   <input_pheno>   <out_dir>   <out_prefix>   [soft_power]");
  q(); 
}
fn_expr      <- as.character(argvs[1])
fn_pheno     <- as.character(argvs[2])
odir         <- as.character(argvs[3])
opref        <- as.character(argvs[4])
if ( is.na(argvs[5]) ) {
	sft_power <- NULL ; 
} else {
	sft_power <- as.character(argvs[5])
}
v_corType <- 'pearson' ; # Could be 'pearson' or 'bicor'; 
# v_corType <- 'bicor'
v_maxPOutliers <- ifelse( v_corType == 'bicor', 0.05, 1 )
# v_robustY      <- ifelse( v_corType == 'bicor', FALSE, TRUE ); # 
v_robustY      <- TRUE ; # Should be FALSE only if v_corType == 'bicor' and trait is binary. 
v_corFnc       <- ifelse( v_corType == 'bicor', 'bicor', 'cor' )

# Fixed parameters; 
moduleTraitPval_thres <- 0.01    # 
v_RsquaredCut     <- c(0.9, 0.8) ; 
v_maxBlockSize    <- 50000 ; # Should be bigger than gene number; 
v_minModuleSize   <- 30 ;    # Minimum gene number required in a module. 
v_mergeCutHeight  <- 0.25 ;  # Cutoff for merging modules. The lower this value is, the fewer modules we get. 
v_networkType     <- 'unsigned'
if ( v_networkType == 'unsigned' || v_networkType == 'signed hybrid' ) {
	powers <- c(1:15);
} else {
	powers <- c(1:30); 
}

### Start working; 
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads() 


ifelse( dir.exists(odir), "exists", dir.create(odir) )
ifelse( file.exists(fn_expr),  "fn_expr OK",  stop(": fn_expr file ", fn_pheno, " not found.\n") )
ifelse( file.exists(fn_pheno), "fn_pheno OK", stop(": fn_pheno file ", fn_pheno, " not found.\n") )
oopref  <- paste0(odir,'/',opref)
fn_sftThres <- paste0(oopref, "_sftThreshold.pdf")
fn_sample_hclust <- paste0(oopref, "_sample_hclust.pdf")
fn_kme  <- paste0(oopref, "_KME.txt")
fn_dendroColorsPng <- paste0(oopref,"_dendroColors.png")
fn_dendroColorsPdf <- paste0(oopref,"_dendroColors.pdf")
fn_adjHmapPng    <- paste0(oopref, "_adjHeatmap.png")
fn_adjHmapPdf    <- paste0(oopref, "_adjHeatmap.pdf")
dir_mebarplot      <- paste0(oopref, "_ME_barplot")
dir_cytoscape  <- paste0(oopref, "_cytoscape/")
fn_ModTrait_HmapPdf <- paste0(oopref, "_module_trait.heatmap.pdf", sep="")
fn_Rimage           <- paste0(oopref, "_hasNet.RData")




RNAseq_expr <- read.table( fn_expr, header=T, row.names=1, sep="\t", na.strings="-" )
# RNAseq_voom <- log2(RNAseq_expr+0.01)
RNAseq_voom <- RNAseq_expr

datExprOri <- t(RNAseq_voom)
( nGenes <- ncol(datExprOri) )
NumberMissingBySample <- apply(is.na(data.frame(datExprOri)),1, sum)
KeepSample <- NumberMissingBySample<0.1*nGenes
table(KeepSample)

datExpr <- datExprOri[KeepSample,]
SampleName <- rownames(datExprOri)[KeepSample]

#SampleName.shrt <- gsub( "Pub\\d_", "", gsub( "self_", "", SampleName ) )
#SampleName.shrt <- gsub("_2012_PI_flesh_D", "_PIFF", SampleName.shrt)
#SampleName.shrt <- gsub("_2012_97_flesh_D", "_97FF", SampleName.shrt)
#SampleName.shrt <- gsub("_2012_97_rind_D", "_97FR",  SampleName.shrt)

nSamples          <- nrow(datExpr)
variancedatExpr   <- as.vector(apply(as.matrix(datExpr),2,var, na.rm=T))
no.missingdatExpr <- as.vector(apply(is.na(as.matrix(datExpr)),2, sum) )
goodExpdatExpr     <- apply( as.matrix(datExpr),2,max ) >= 1
KeepGenes         <-  variancedatExpr>0 & no.missingdatExpr<0.1*nSamples & goodExpdatExpr
table(KeepGenes)

datExpr    <- datExpr[, KeepGenes]
GeneName   <- colnames(datExpr)
save(datExpr, file=paste0(oopref, "_datExprInput.RData"))


# tree       <- hclust(dist(datExpr),method = 'average')
tree       <- hclust(dist(log2(datExpr+0.01)),method = 'average') # The log2(exp+0.01) produces more reasonable result. 
pdf(file= fn_sample_hclust, , width= 7, height= 5)
plot(tree,xlab="", sub="", cex = 0.7,main="Sample Clustering")
plot(tree,xlab="", sub="", cex = 0.7,main="Sample Clustering", hang=-1)
dev.off()

if ( is.null(sft_power) ) {
  # https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
  sft = pickSoftThreshold(
    datExpr, 
    RsquaredCut = v_RsquaredCut[1],
    powerVector = powers, 
    networkType = v_networkType, 
    verbose = 5
  )
  
  if ( is.na(sft$powerEstimate) ) {
    # http://blog.genesino.com/2018/04/wgcna/#经验power-无满足条件的power时选用
    if ( length(v_RsquaredCut) > 1 ) {
      test_RsquaredCut <- v_RsquaredCut[-1]
      for ( v1 in test_RsquaredCut ) {
        for ( v2 in 1:length(sft$fitIndices$Power) ) {
          if ( sft$fitIndices$SFT.R.sq[v2] >= v1 ) {
            message("Soft power determined by threshold: ", sft$fitIndices$SFT.R.sq[v2], " as ", sft$fitIndices$Power[v2]); 
            sft_power <- sft$fitIndices$Power[v2]; 
            break; 
          }
	}
        if ( !is.null(sft_power) ) {
          break; 
        }
      }
    }
    if ( is.null(sft_power) ) {
      if ( nSamples < 20 ) {
        sft_power <- ifelse( v_networkType == "unsigned" || v_networkType == "signed hybrid", 9, 18 ); 
        message("Soft power determined by empirical value with [nSamples, networkType] [", nSamples, ",", v_networkType, "] as : ", sft_power); 
      } else if ( nSamples < 30 ) {
        sft_power <- ifelse( v_networkType == "unsigned" || v_networkType == "signed hybrid", 8, 16 ); 
        message("Soft power determined by empirical value with [nSamples, networkType] [", nSamples, ",", v_networkType, "] as : ", sft_power); 
      } else if ( nSamples < 40 ) {
        sft_power <- ifelse( v_networkType == "unsigned" || v_networkType == "signed hybrid", 7, 14 ); 
        message("Soft power determined by empirical value with [nSamples, networkType] [", nSamples, ",", v_networkType, "] as : ", sft_power); 
      } else {
        sft_power <- ifelse( v_networkType == "unsigned" || v_networkType == "signed hybrid", 6, 12 ); 
        message("Soft power determined by empirical value with [nSamples, networkType] [", nSamples, ",", v_networkType, "] as : ", sft_power); 
      }
    }
  } else {
    sft_power <- sft$powerEstimate
  }

  pdf( file=fn_sftThres, width=7, height=4 )
  # sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  plot(
    sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)",
    ylab="Scale Free Topology Model Fit,signed R^2",
    type="n", main = paste("Scale independence"),
    ylim=c(-1,1)
  );
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red"); 
  abline(h=0.90,col="red")
  abline(h=0.80,col="blue")
  legend("bottomright", legend=c(paste0("threshold0=", sft$powerEstimate), paste0('threshold1=', sft_power)))
  plot(
    sft$fitIndices[,1], sft$fitIndices[,5], 
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity")
  ); 
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  # dev.copy2pdf(file=fn_sftThres, width=10, height=10)
  dev.off()

}

#### 
message("Soft power threshold used is :", sft_power)
net = blockwiseModules(
  datExpr,
  power = sft_power,                          #软阈值，前面计算出来的, beta值
  maxBlockSize = v_maxBlockSize,              #最大block大小，应>=基因总数，将所有基因放在一个block中, 涉及内存使用. 
  TOMType = "unsigned",                       #选择unsigned，使用标准TOM矩阵
  networkType = v_networkType,                #signed recommended. http://www.peterlangfelder.com/signed-or-unsigned-which-network-type-is-preferable/
  deepSplit = 2,                              #剪切树参数，deepSplit取值0-4
  minModuleSize = v_minModuleSize,            # 通常最小值设置为30, 但也可以试试10/20之类的, 每次模块的内容都会发生变化; 
  mergeCutHeight = v_mergeCutHeight,          # 模块合并参数，越大模块越少
  numericLabels = TRUE,                       # T返回数字，F返回颜色
  pamRespectsDendro = FALSE,  
  saveTOMs = TRUE,                            # 保存TOM矩阵(放到当前目录), 以便以后重用
  saveTOMFileBase = paste0(oopref, "_TOM"),   # 保存TOM矩阵的文件的前缀
  loadTOMs = TRUE,                            # 导入生成的TOM矩阵; 
  corType      = v_corType,                   # 'bicor' recommended. 
  maxPOutliers = v_maxPOutliers, 
  verbose = 3
)
save.image( file= fn_Rimage )

moduleColors    <- labels2colors(net$colors)
moduleLabels    <- net$colors
MEsList         <- moduleEigengenes(datExpr, moduleColors) ; 
MEs             <- orderMEs(MEsList$eigengenes); 
rownames( MEs ) <- SampleName 
net.MEs         <- net$MEs ;
net.geneTree    <- net$dendrograms[[1]] ;
save( moduleColors, moduleLabels, MEs, net.MEs, net.geneTree, file= paste0(oopref, "_network_ForPhenoAssoc.RData") )

# if ( file.exists(fn_kme) ) { file.remove(fn_kme) }
t.KME           <- signedKME( datExpr, MEs, corFnc=v_corFnc, corOptions = "use = 'p'", outputColumnName= "ME" ); # This is genes' module membership (MM); 
t.ModuleID      <- moduleColors
t.ModuleIDME    <- paste0("ME", t.ModuleID)
t.v1 <- data.frame( MC=match(t.ModuleIDME, colnames(t.KME)), t.KME ) ; 
t.KMEinM    <- apply(t.v1, 1, FUN=function(x){ x[x[1]+1] }); 
t.KMEinMabs <- abs( t.KMEinM )
colnames(t.KME) <- paste0("gMM_", substring(colnames(t.KME), 3))
t.connect       <- intramodularConnectivity.fromExpr(datExpr, moduleColors, corFnc=v_corFnc, corOptions = "use = 'p'", networkType=v_networkType )
kme_tbl         <- data.frame(EleID=colnames(datExpr), KME=t.KMEinM, ModuleID=t.ModuleID, absKME=t.KMEinMabs, ModuleNum=net$colors, t.connect[,c("kWithin","kOut")], t.KME)
kme_tbl         <- kme_tbl[order(kme_tbl$ModuleNum, -kme_tbl$absKME),]
write.table( kme_tbl, file=fn_kme, col.names=T, row.names=F, quote=F, sep="\t", append= F )
rm( list=strsplit("t.KME t.ModuleID t.ModuleIDME t.v1 t.KMEinM t.KMEinMabs", " ")[[1]] )


expColor           <- t(numbers2colors(log2(datExpr+0.01),colors=blueWhiteRed(100),naColor="grey"))
colnames(expColor) <- rownames(datExpr)
# png( fn_dendroColorsPng ,height = 700,width = 900 )
pdf( fn_dendroColorsPdf , height= 14, width= 14 )
plotDendroAndColors(
  net$dendrograms[[1]], 
  colors=cbind(moduleColors[net$blockGenes[[1]]],expColor),
  c("Module",colnames(expColor)),
  dendroLabels = F, hang = 0.03,
  addGuide = T, guideHang = 0.05,
  cex.rowText=0.5
)
dev.off()

# png( fn_adjHmapPng,height = 1000,width = 900 )
pdf( fn_adjHmapPdf, height = 7,width = 7 )
plotEigengeneNetworks(
  multiME    = MEs, 
  setLabels  = "Eigengene adjacency", 
  plotDendrograms = T, 
  marDendro  = c(4,4,2,4), 
  marHeatmap = c(4,10,2,10)
)
dev.off()

# 
ifelse( dir.exists(dir_mebarplot), "exists", dir.create(dir_mebarplot) )
for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  ME=MEs[,paste("ME",module,sep="")]
  png(paste0(dir_mebarplot, "/", module, ".express.barplot.png", sep=""), height = 700,width = 900)
  # pdf(paste0(dir_mebarplot, "/", module, ".express.barplot.pdf", sep=""), height= 7, width= 7 )
  par(mfrow=c(2,1),mar=c(0.3,6.5,5,3))
  plotMat(t(scale(datExpr[,moduleColors==module])),
  # rlabels=F,main=module,cex.main=2,clabels=F)
  rlabels=F,main=module,cex.main=2,clabels=SampleName)
  
  par(mar=c(5,4.2,0,0.7))
  barplot(ME,col=module,main="",cex.main=2,ylab="eigengene expression",xlab="sample")
  dev.off()
}

# For cytoscape : 
ifelse( dir.exists(dir_cytoscape), "exists", dir.create(dir_cytoscape) )
### wgcna_dat1/dat1_TOM-block.1.RData 
load( paste0(oopref, "_TOM-block.1.RData") ) ; # This is only required if we are redoing this work. 
#TOM = TOMsimilarityFromExpr(datExpr, power =sft_power,TOMType = "unsigned", networkType=v_networkType); 
TOM.mat <- as.matrix(TOM)

for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, module))
  modProbes = probes[inModule]
  modTOM = TOM.mat[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  cyt = exportNetworkToCytoscape(
          modTOM,
          # edgeFile = paste("CytoscapeInput-edges-", module, ".txt", sep=""),
          edgeFile = paste0(dir_cytoscape, "input-edges-", module, ".txt", sep=""), 
          # nodeFile = paste("CytoscapeInput-nodes-", module, ".txt", sep=""),
          nodeFile = paste0(dir_cytoscape, "input-nodes-", module, ".txt", sep=""), 
          weighted = TRUE,
          threshold = 0.02,
          nodeNames = modProbes,
          nodeAttr = moduleColors[inModule]
        )
}

# Association : Here I want to add each sample as a pheno. 
### 'NA' in fn_pheno is accepted with WGCNA::cor(use="pairwise.complete.obs"), but this tends to provide lower p-values than using zero. 
### Varaibles required from previous analysis : 
##### v_corType
##### fn_pheno 
##### fn_ModTrait_HmapPdf
##### moduleTraitPval_thres
##### MEs     : Get from paste0(oopref, "_network_ForPhenoAssoc.RData"); 
##### SampleName # Should be rownames(MEs) or rownames(datExpr)
##### nSamples   # Should be nrow(MEs) or nrow(datExpr)
### So I want to save data for association anaylysis. 

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
    corP.mod_self              <- WGCNA::bicor(MEs_Ass, use = "pairwise.complete.obs")
    cor.mod_self               <- corP.mod_self$bicor
    pval.mod_self              <- corP.mod_self$p
    corP.mod_trait             <- WGCNA::bicor(MEs_Ass, datTraits, use = "pairwise.complete.obs", robustY = v_robustY)
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

save.image( file= fn_Rimage )

message("All done"); 

