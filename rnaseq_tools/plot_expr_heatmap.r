#!/home/Sunhh/bin/Rscript
argvs <- commandArgs( trailingOnly=TRUE );
if ( is.na(argvs[1]) ) {
  message("Rscript this.R   <exp_mat.tab>  <num_of_group_to_divide>\nOutput: exp_mat.tab.pdf, exp_mat.tab.grp\n");
  q();
}

fn_expMat        <- as.character(argvs[1]) ;
numGrp           <- as.integer(argvs[2]) ;

fn_outPdf <- paste0(fn_expMat, ".pdf")
fn_outGrp <- paste0(fn_expMat, ".grp")

# [No need to change]
#   Load the expression matrix.
expMat <- read.table(fn_expMat, header=T, row.names=1, stringsAsFactors=F)
#   Use log(1+expression) values instead of the raw RPKM values for clustering.
logExp <- log1p(expMat)

# [No need to change] Load the libraries required.
library(RColorBrewer)
library(pheatmap)

# Set the colors used in the heatmap plotting.
# "RdBu" can be replaced by: Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd BrBG PiYG PRGn PuOr RdGy RdYlBu RdYlGn Spectral
# No need to change other values.
cc <- colorRampPalette( rev(brewer.pal(n= 9, name= "RdBu")) )(255)

# Plot the heatmap figure.
pdf(file=fn_outPdf, width=8, height=8)
h=pheatmap(logExp,
  main = "Expression heatmap", # Set the title of the plot.
  scale= 'row',     # Further center the expression matrix within each gene across samples.
  na_col = "grey", color=cc, # Assign the colors used in the heatmap.
  cluster_rows =T,  # Value "T" means to cluster the rows (genes), and "F" will turn off this clustering.
  show_rownames = F,# "F" means not to show the name of rows (genes) on the plot.
  cluster_cols = T, # Value "T" means to cluster the columns (samples), and "F" turns off this.
  show_colnames=T,  # "T" means to show the name of columns (samples) on the plot.
  display_numbers=F,# "F" means not to plot the expression values in the plot.
  legend=T         # "legend=T" means to show the scale bar of scaled log(1+expression) values on the right.
)
dev.off()

# Write the grouping information.
memb <- cutree(h$tree_row, k = numGrp)
memb1 <- data.frame(memb)
memb1$Gene_ID <- names(memb)
memb1 <- memb1[,c(2,1)]
write.table(memb1, file=fn_outGrp, sep="\t", quote=F, col.names= c("Gene_ID", "Group_ID"), row.names=F )


