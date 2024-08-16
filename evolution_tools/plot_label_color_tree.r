library(ggtree)
library(dplyr)
library(ggplot2)
library(tidytree)

##############################
# Load the ID mapping file.
##############################
### Change oldID to newID on the tree.
fn_IDmap  <- 'GSID_oldID_newID_groupID'
dta_IDmap <- dplyr::as_tibble ( read.table( fn_IDmap, header=T, stringsAsFactor=F ) ) %>% dplyr::mutate( label= oldID ) %>% dplyr::filter( !is.na(label) )


##############################
# Load the newick-format tree.
##############################
fn_nwk <- 'in.nwk'
tree0 <- treeio::read.newick( fn_nwk, node.label= 'support' ) # Store bootstrap value as 'support'.
tree0_tbl <- as_data_frame( tree0 ) %>% dplyr::mutate( bts= round( support * 100, 0) ) # Multiply bootstrap value with 100 as 'bts'.

##############################
# Set up output file names.
##############################
fo_treePdf <- 'out.pdf'

##############################
# Modify the tree as needed.
##############################

### Re-root tree.
tree1 <- tree0
##### With a single individual 'outgroup_INDV' node as outgroup.
# tree1 <- ggtree::reroot( as.phylo(tree0), node= nodeid( tree0, 'outgroup_INDV') )
##### With the parent node of 'outgroup_INDV1' as the outgroup.
# tree1 <- ggtree::reroot( as.phylo(tree0), node= phytools::getParent( as.phylo(tree0), node= nodeid( tree0, 'outgroup_INDV1') ) )

### Rename IDs on the tree.
##### Map tree IDs with node IDs.
tree1_tbl <- as_data_frame( tree1 )
##### Set 'node_map' attribute to tree1 variable (attr()), then convert the table into tbl_df object (as_tibble()).
tree1_o2n <- attr( tree1, 'node_map' ) %>% dplyr::as_tibble()
##### Add tree1_o2n table to tree1_tbl according to node mapping ('to').
tree1_tbl_j1 <- dplyr::mutate( tree1_tbl, to= node ) %>% left_join( y= tree1_o2n, by= 'to' )
##### Further include tree0_tbl.
tree1_tbl_j2 <- left_join( x= tree1_tbl_j1, y= dplyr::mutate( tree0_tbl, from=node ), by= 'from' ) %>% dplyr::mutate( label= label.x )
##### Further include newIDs (label) from dta_IDmap.
jn_ID <- left_join( x=tree1_tbl_j2, y=dta_IDmap, by= 'label' )
##### Replace the IDs in tree1.
jn_ID_tip <- jn_ID[ !is.na(jn_ID$label), ]
tree1$tip.label <- as.character( jn_ID_tip$newID )

### Set colors for taxa groups.
##### Use 'groupID' column to divide groups.
grp_cc <- as.factor( jn_ID_tip$groupID )
##### Assign colors for each group.
basic_colors <- c(
  rgb( maxColorValue=255, 0, 100, 255), 
  rgb( maxColorValue=255, 255, 0, 0), 
  rgb( maxColorValue=255, 0, 200, 0), 
  rgb( maxColorValue=255, 200, 0, 255), 
  rgb( maxColorValue=255, 160, 160, 160), 
  rgb( maxColorValue=255, 255, 160, 0), 
  rgb( maxColorValue=255, 0, 0, 0),
  rgb( maxColorValue=255, 0,100,0)
)
levels( grp_cc ) <- basic_colors[1:length(levels(grp_cc))]


##############################
# Plot the final tree.
##############################

p0 <- ggtree ( tree1 , layout= 'fan', ladderize=FALSE, branch.length='none' )
p1 <- p0 + geom_tippoint( color= as.character( grp_cc )  ) 
# p1 + geom_text2( aes( subset=!isTip, label=c(jn_ID$bts)) )
p2 <- p1 + geom_nodelab( aes( label= c(jn_ID$bts) ), geom= 'text', hjust=0.5 ) + geom_tiplab2()

### Rotate tree.
p2 <- rotate_tree(p2, angle=90)

### Flip two clades.
child_node_1 <- 'child_1_INDV_node'
child_node_2 <- 'child_2_INDV_node'
mrca_node <- ggtree::MRCA( tree1, tip= c( child_node_1, child_node_2))
mrca_child_nodes <- phytools::getDescendants( as.phylo(tree1), node= mrca_node )
p2 <- flip( p2, mrca_child_nodes[1], mrca_child_nodes[2] ) %>% rotate_tree( angle=90 )

pdf(file= fo_treePdf, width=21, height=21)
print(p2)
dev.off()



