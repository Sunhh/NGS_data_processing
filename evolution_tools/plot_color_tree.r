library(ggtree)
library(cowplot)
library(ggplot2)  
library(dplyr)
library(ggnewscale)
library(tidyverse)
library(ape)
library(data.table)
library(RColorBrewer)

 
# prune
tree <- read.tree("emQ30-4dA-treID.nwk")
write.table(tree$tip.label, "emQ30-4dA-treID.nwk.label.txt", quote=FALSE, sep="\t", col.name=FALSE, row.name=FALSE)

# assign origin color according to order in tree$tip.label
# read the group information
group <- read.table("emQ30-4dA-treID.nwk.label.txt.pop1", header = FALSE)

my_info <- data.table(tip_lbs = tree$tip.label, groupC = group[,2])

grC <- split(my_info$tip_lbs, my_info$groupC)

 

# use groupOTU for taxa classification

tree_grC <- ggtree::groupOTU(tree, grC)

 

# Define the group color

my_cols <- c("#FF0000","#548235",
             "#C65911","#7030A0",
             "#DFE349","#00B0F0",
             "#D400E3","#FFC000", "#000000")

my_cols <- my_cols[1:length(unique(my_info$groupC))]
names(my_cols) <- levels(attributes(tree_grC)$group)

# Check the group color

#scales::show_col(my_cols); my_cols

 

pdf("emQ30-4dA-treID-col_1.pdf", height = 32, width = 45)

# circular tree
ggtree(tree_grC, layout="circular", branch.length="none",mapping = aes(color = group))+

# rectangular tree
#ggtree(tree_grC, branch.length="none", mapping = aes(color = group))+  
geom_tiplab(size = 3, aes(angle=angle)) +

#geom_tiplab(size = 3) + scale_y_reverse() + 
scale_color_manual(name = 'Taxonomy', values = my_cols) +

theme(legend.title=element_text(size=15), legend.text=element_text(size=12))

dev.off()


