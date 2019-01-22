#!/home/Sunhh/bin/Rscript
# 2018-11-05 

# options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# source("http://bioconductor.org/biocLite.R")
# biocLite("topGO")
# biocLite('Rgraphviz')

# Setup basic options. 
fdr_thres      <- 0.01 ; 
v_nodeSize     <- 3 ; 
v_toShowN <- 10

argvs <- commandArgs( trailingOnly=TRUE ); 
if ( is.na(argvs[1]) ) { 
  message("Rscript this.R   <sepTbl.genID_GOIDorNA>   <DEG.genID_FDR>   <out_prefix>   [special_GO.list]");
  message("\n\nOutput : \n  opref.GO_BP/MF/CC.txt/pdf\n"); 
  q();
}

fn_geneID2GO   <- as.character( argvs[1] ); 
# fn_geneID2GO   <- '/Data/Sunhh/database/db_fasta/watermelon/97103/pb/v2/Final/annot/for_GO/ordered_feiID.gene2GO_topGO_20181105.sepTbl'
fn_geneSubset  <- as.character( argvs[2] ); 
# fn_geneSubset  <- 'DEG_dat1_in_wwW_FF'  # Should be three column : GeneID \\t FDR_NA \\t DEG_txt_'H/L/N'
opref          <- as.character( argvs[3] ); 
# opref          <- fn_geneSubset
if ( is.na(argvs[4]) ) {
  fn_specGO    <- NULL; 
  # fn_specGO      <- '/Data/Sunhh/database/db_fasta/watermelon/97103/pb/v2/Final/annot/for_GO/specialClass/ref/go_class.txt'
} else {
  fn_specGO    <- as.character( argvs[4] )
}

### Codes; 
library(topGO)
library(Rgraphviz)
deg_fdr <- function( fdrScore ) { return( fdrScore < fdr_thres )  }

# Load in GO annotation and all gene names. 
### geneID2GO <- readMappings( file= fn_geneID2GO ); # Here 
geneID2GO_all <- read.table( fn_geneID2GO, header=F, stringsAsFactors=F, sep="\t", quote="" )
k.hasGO   <- !is.na( geneID2GO_all[,2] )
geneID2GO <- by( geneID2GO_all[k.hasGO,2], geneID2GO_all[k.hasGO,1], function(x) as.character(x) )
goID2gene <- inverseList( geneID2GO )
geneNames_all <- unique( geneID2GO_all[,1] )

# Load in gene subset list. 
data_subset  <-read.table(fn_geneSubset, row.names = NULL, header=FALSE, check.names =F, stringsAsFactors=F, quote="" )
if ( ncol(data_subset) == 1 ) { data_subset[,2] <- 0 } # Add geneScore (FDR)
if ( ncol(data_subset) == 2 ) { data_subset[,3] <- 'H' } # Add DEG_txt ; 
k.na <- is.na( data_subset[,2] )
data_subset[k.na, 2] <- 1  
data_subset[k.na, 3] <- 'N'
if ( nrow(data_subset) == length(geneNames_all) ) {
  geneScore <- data_subset[,2] ; # +1 ?
  names(geneScore) <- data_subset[,1]
} else {
  geneScore <- factor( as.integer(geneNames_all %in% data_subset[,1]) )
  names(geneScore) <- geneNames_all
}

# Load in special GO terms; 
if (!is.null(fn_specGO)) {
  data_specGO <- read.table( fn_specGO, header=F, sep="\t", stringsAsFactors=F, quote="" )
  if ( data_specGO[1,1] == "GO_ID") {
    data_specGO <- data_specGO[-1, ]
  }
  specGOs <- data_specGO[,1]
}

# GO enrichment; 
# https://www.jianshu.com/p/9e21f2196178
# https://www.cnblogs.com/djx571/p/9625830.html
# The output are raw p-value instead of FDR !!!
for ( topTerm in c('BP', 'MF', 'CC') ) {
  sampleGOdata <- new(
    "topGOdata", 
    ontology   = topTerm, 
    allGenes   = geneScore, 
    annot      = annFUN.gene2GO, 
    gene2GO    = geneID2GO, 
    nodeSize   = v_nodeSize, 
    geneSelectionFun = deg_fdr
  )
  # whichTests() ; whichAlgorithms() ; 
  resultFis.cla    <- runTest( sampleGOdata, algorithm= 'classic', statistic = 'fisher' ) # Fisher’s exact test, based on gene counts; 
  resultKS.elim    <- runTest( sampleGOdata, algorithm= 'elim',    statistic = 'ks' ) ; # Kolmogorov-Smirnov like test, based on gene scores; 
  resultKS.cla     <- runTest( sampleGOdata, algorithm= 'classic', statistic = 'ks' ) ; # Kolmogorov-Smirnov like test, based on gene scores; 
  nodeN            <- length( score( resultFis.cla ) ) ; 
  allRes           <- GenTable( 
                         sampleGOdata, 
                         classicFisher = resultFis.cla, 
                         elimKS        = resultKS.elim, 
                         classicKS     = resultKS.cla, 
                         orderBy       = 'classicFisher', 
                         ranksOf       = 'classicFisher', 
                         topNodes      = nodeN
  )
  toShowN <- v_toShowN
  if ( nodeN < toShowN ) { toShowN <- nodeN }
  v1 <- apply(as.matrix(allRes$GO.ID), MARGIN=1, FUN=function(x) { paste0(unlist(genesInTerm(sampleGOdata, whichGO=x[1])), collapse=",") } )
  write.table( cbind(allRes, Gene.IDs=v1), file= paste0(opref, ".GO_", topTerm, ".txt", sep=''), sep="\t", quote=F, col.names=T, row.names=F )
  pdf( file=paste0(opref, ".GO_", topTerm, ".topNodes.pdf"), width=10, height=10 )
  showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = toShowN, useInfo = "all")
  legend('topleft', legend=c( paste0(opref, ".GO_", topTerm, ".elimKS.top", toShowN, "nodes" ) ) )
  showSigOfNodes(sampleGOdata, score(resultFis.cla), firstSigNodes = toShowN, useInfo = "all")
  legend('topleft', legend=c( paste0(opref, ".GO_", topTerm, ".classicFisher.top", toShowN, "nodes") ) )
  showSigOfNodes(sampleGOdata, score(resultKS.cla), firstSigNodes = toShowN, useInfo = "all")
  legend('topleft', legend=c( paste0(opref, ".GO_", topTerm, ".classicKS.top", toShowN, "nodes") ) )
  dev.off() 
  
  if (!is.null(fn_specGO)) {
    specGOs_has <- sampleGOdata@graph@nodes[ sampleGOdata@graph@nodes %in% specGOs ]
    if ( length(specGOs_has) > 0 ) {
      pdf( file=paste0(opref, ".GO_", topTerm, ".specGOs.pdf") , width=10, height=10)
      showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = NULL, wantedNodes = specGOs_has, useInfo = "all")
      legend('topleft', legend=c( paste0(opref, ".GO_", topTerm, ".elimKS.specGOs") ) )
      showSigOfNodes(sampleGOdata, score(resultFis.cla), firstSigNodes = NULL, wantedNodes = specGOs_has, useInfo = "all")
      legend('topleft', legend=c( paste0(opref, ".GO_", topTerm, ".classicFisher.specGOs") ) )
      showSigOfNodes(sampleGOdata, score(resultKS.cla), firstSigNodes = NULL, wantedNodes = specGOs_has, useInfo = "all")
      legend('topleft', legend=c( paste0(opref, ".GO_", topTerm, ".classicKS.specGOs") ) )
      dev.off() 
    }
  }
}
# bp.allRes.KS_elim   <- GenTable( bp.sampleGOdata, KS= bp.resultKS.elim,      ranksOf= 'classic', topNodes= attributes( bp.resultKS.elim )$geneData[4] )
# bp.allRes.Fis_cla   <- GenTable( bp.sampleGOdata, classic= bp.resultFis.cla, ranksOf= 'classic', topNodes= attributes( bp.resultFis.cla )$geneData[4] )
# pValue.KS.classic <- score(bp.resultKS.cla)
# pValue.KS.elim <- score(bp.resultKS.elim)[names(pValue.KS.classic)]
# gstat <- termStat(bp.sampleGOdata, names(pValue.KS.classic))
# gSize <- gstat$Annotated / max(gstat$Annotated) * 4
# plot(pValue.KS.classic, pValue.KS.elim, xlab = "p-value classic", ylab = "p-value elim",pch = 19, cex = gSize)

# BP part object: 
# genes( bp.sampleGOdata ) ; # Get annotated genes; 
# description( bp.sampleGOdata ) ; # Get description of this data; 
# geneScore( bp.sampleGOdata, whichGenes= names(geneList), use.names=FALSE ) ; # Get score of gene list; 
# graph( bp.sampleGOdata ) ; # Get graph information; 
# usedGO( bp.sampleGOdata ) ; # Get all used GO terms; 
# termStat( bp.sampleGOdata, selected.GOterms ) ; # According to selected GO term list (selected.GOterms) to stat genes; 
# ggg <- integer(length(geneNames_all))+1
# ggg[ geneNames_all %in% genes( bp.sampleGOdata ) ] <- 0
# ggg[ !( geneNames_all %in% data_subset[,1] ) ] <- 2
# ggg <- factor( ggg, labels=c('Used', 'Not annotated', 'Filtered') )
# table(ggg)
# pValues <- getPvalues(exprs(ALL), classlabel=y, alternative="two.sided")
# geneVar <- apply(exprs(ALL), 1, var)
# dd <- data.frame(x = geneVar[allGenes], y=log10(pValues[allGenes]), groups=group)
# lattice::xyplot(y~x|group, data=dd, groups=group)


message("All done."); 

