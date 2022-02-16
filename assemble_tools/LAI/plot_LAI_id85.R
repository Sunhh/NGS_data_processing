# [2/16/2022] Need to modify the code to fit 'ChrS' value at each run.

tblF   <-   'pepper_genome.fa.out.LAI'
chrF   <-   'pepper_genome.fa.chrID'
# gapF   <-   'WM97pbV1_fwd.chr.fa.n_list'
outPdf <-   'pepper_genome.fa.out.LAI_unlock_distr_byChr.pdf'

vv     <-   'LAI'

library(dplyr)
library(ggplot2)

# Input data : 
aa    <- read.table( tblF, header=T, stringsAsFactors=F, sep="\t" )
cc    <- read.table( chrF, header=F, stringsAsFactors=F, sep="\t" )
# nlist <- read.table( gapF, header=T, stringsAsFactors=F, sep="\t" )

ds <- aa %>% tibble::as_tibble() %>% dplyr::filter( Chr %in% cc[,1] )


# Prepare data for plotting; 'vUse' is the data to be plotted; 
t1 <- factor( ds$Chr )
levels(t1) <- rep_len(c("blue", "red"), length.out=length(levels(t1)))
toPlot <- dplyr::select( ds, one_of('Chr', 'From', eval(vv)) ) %>% dplyr::mutate( coll = as.character(t1), ChrS= substr(Chr, 4, 5) )
colnames(toPlot)[ colnames(toPlot) %in% vv ] <- 'vUse'



# Plot 

pdf( file= outPdf, width=14, height=3.5 ) 
p <- ggplot2::ggplot( data= toPlot )
p <- p + ggplot2::facet_grid(. ~ ChrS , scales= 'free_x' , space = 'free_x' )
p <- p + ggplot2::labs( y= vv )
p1 <- p + theme( panel.grid =element_blank(), panel.background=element_blank(), axis.line= element_line(colour= "black"), axis.ticks.x= element_blank(), axis.text.x= element_blank(), axis.title.x= element_blank() )
# p <- p + ggplot2::geom_point( mapping= ggplot2::aes(x= From, y= vUse, colour= colour ) )
p2 <- p1 + ggplot2::geom_point( mapping= ggplot2::aes(x= From, y= vUse), colour= toPlot$coll, size=0.2 )
p2
dev.off()


