#!/home/Sunhh/bin/Rscript
argvs <- commandArgs( trailingOnly=TRUE );
if ( is.na(argvs[1]) ) {
  message("Rscript this.R   <sim_mat.tab>   <out_pdf_file>");
  q();
}
fn_simMat        <- as.character(argvs[1]) ; # Out put of : perl dist_of_twoKME.pl sign redo/wgcna_dat1_signed/dat1_KME.txt unsi redo/wgcna_dat1_unsign/dat1_KME.txt > tt
fn_outPdf        <- as.character(argvs[2]) ; 

library(pheatmap)
library(RColorBrewer)
colors <- colorRampPalette( rev(brewer.pal(n= 9, name= "OrRd")) )(255)

simMat <- read.table(fn_simMat, header=T, row.names=1, stringsAsFactors=F)
distMat <- as.dist(1-simMat)
distMat.mat <- as.matrix(distMat)

pdf(file=fn_outPdf, width=8, height=8)
pheatmap(distMat.mat, 
  clustering_distance_rows = distMat , 
  clustering_distance_cols = distMat , 
  col = colors , 
  show_colnames=TRUE)
dev.off()

