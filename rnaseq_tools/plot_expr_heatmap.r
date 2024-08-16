#!/home/Sunhh/bin/Rscript

##############################
# Setup working directory.
##############################
# setwd("D:\\Work\\Watermelon\\traits\\gynoecious\\plot")

##############################
# Load libraries.
##############################
# [No need to change] Load the libraries required.
# install.packages('devtools')
#library(devtools)
#install_github("raivokolde/pheatmap")
library(RColorBrewer)
library(pheatmap)

##############################
# Load input files and parameters: input matrix file, number of groups to cluster the genes.
##############################
argvs <- commandArgs( trailingOnly=TRUE );
if ( is.na(argvs[1]) ) {
  message("\nRscript this.R   <exp_mat.tab>  <num_of_group_to_divide>\nOutput: exp_mat.tab.pdf, exp_mat.tab.grp\n");
  message("##############################
# Example of input file <exp_mat.tab>
##############################
# ID           GG_XHB   gg_XHBGM  GG_ZZJM  gg_T1
# Cla001867    49.16    103.02    5.05     34.47
# Cla002628    10.35    22.8      4.06     10.64
# Cla002678    10.44    21.67     7.4      28.42
# Cla003357     3.37     7.42     5.49     21.8
# Cla003741    14.82    30.29     9.35     42.59
");
  q();
}

fn_expMat        <- as.character(argvs[1]) ;
numGrp           <- as.integer(argvs[2]) ;

fn_outPdf <- paste0(fn_expMat, ".pdf")
fn_outGrp <- paste0(fn_expMat, ".grp")

# [No need to change]
#   Load the expression matrix.
expMat <- read.table(fn_expMat, header=T, row.names=1, stringsAsFactors=F)
# expMat <- read.delim(header=T,file= fn_expMat, row.names = 1,sep="\t")
#   Use log(1+expression) values instead of the raw RPKM values for clustering.
logExp <- log1p(expMat)

##############################
# Set the colors used in the heatmap plotting.
# "RdBu" can be replaced by: Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd BrBG PiYG PRGn PuOr RdGy RdYlBu RdYlGn Spectral
##############################
# No need to change other values.
cc <- colorRampPalette( rev(brewer.pal(n= 9, name= "RdBu")) )(255)

##############################
# Plot the heatmap figure into files.
##############################
### Missing data in 'logExp' matrix may alter the clustering or even make it failed.
### Use 'scale' to determine how to normalize the data.
pdf(file=fn_outPdf, width=8, height=8)
h=pheatmap(logExp,
  main = "Expression heatmap", # Set the title of the plot.
  scale= 'row',                # Further center/normalize the expression matrix within each gene across samples.
  na_col = "grey", color=cc,   # Assign the colors used for NA and groups in the heatmap.
  cluster_rows =T,             # Value "T" means to cluster the rows (genes), and "F" will turn off this clustering.
  show_rownames = F,           # "F" means not to show the name of rows (genes) on the plot.
  cluster_cols = T,            # Value "T" means to cluster the columns (samples), and "F" turns off this.
  show_colnames=T,             # "T" means to show the name of columns (samples) on the plot.
  # cellwidth= 20, cellheight= 20,
  display_numbers=F,           # "F" means not to plot the expression values in the plot.
  legend=T                     # "legend=T" means to show the scale bar of scaled log(1+expression) values on the right.
)
dev.off()

##############################
# Write the grouping information.
##############################
memb <- cutree(h$tree_row, k = numGrp)
memb1 <- data.frame(memb)
memb1$Gene_ID <- names(memb)
memb1 <- memb1[,c(2,1)]
write.table(memb1, file=fn_outGrp, sep="\t", quote=F, col.names= c("Gene_ID", "Group_ID"), row.names=F )


