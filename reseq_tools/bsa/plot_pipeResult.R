#!/usr/bin/Rscript

argvs <- commandArgs( trailingOnly=TRUE ) ; 
# Example of command : 
#   Rscript plot_Gprime.R set1_red2whiteFlesh.bsaPipe.w100k chrList Gprime set1_red2whiteFlesh.bsaPipe.w100k.Gprime.pdf

tblF <- as.character( argvs[1] ) ;
chrF <- as.character( argvs[2] )
vv   <- as.character( argvs[3] )
outPdf <- as.character( argvs[4] )

library(dplyr)
library(ggplot2)
format_genomic <- function(...) {
  # Based on code by Ben Tupper
  # https://stat.ethz.ch/pipermail/r-help/2012-January/299804.html
  # Args:
  #   ...: Args passed to format()
  #
  # Returns:
  #   A function to format a vector of strings using
  #   SI prefix notation
  #
  function(x) {
    limits <- c(1e0,   1e3, 1e6)
    prefix <- c("","K","M")

    # Vector with array indices according to position in intervals
    i <- findInterval(abs(x), limits)

    # Set prefix to " " for very small values < 1e-24
    i <- ifelse(i==0, which(limits == 1e0), i)

    paste0( format(round(x/limits[i], 1), trim=TRUE, scientific=FALSE, ...), prefix[i] )
  }
}

# Input data : 
aa <- read.table( tblF, header=T, stringsAsFactors=F )
cc <- read.table( chrF, header=F, stringsAsFactors=F )
if ( sum(colnames(aa) == "CHROM") == 1 ) {
	; 
} else if ( colnames(aa)[1] == "ChromID" ) {
	colnames(aa)[1] = "CHROM"
}
if ( sum(colnames(aa) == "POS") == 1 ) {
	; 
} else if ( colnames(aa)[2] == "WindS" ) {
	colnames(aa)[2] = "POS"
}
ds <- aa %>% tibble::as_tibble() %>% dplyr::filter( CHROM %in% cc[,1] )

# Prepare data for plotting; 'vUse' is the data to be plotted; 
toPlot <- dplyr::select( ds, one_of('CHROM', 'POS', eval(vv)) )
colnames(toPlot)[ colnames(toPlot) %in% vv ] <- 'vUse'

# Plot 
pdf( file= outPdf, width=21, height=3.5 )
p <- ggplot2::ggplot( data= toPlot )
p <- p + 
	ggplot2::geom_line( mapping= ggplot2::aes(x= POS, y= vUse) ) +
	ggplot2::labs( y= vv ) +
	ggplot2::facet_grid(
		. ~ CHROM ,
		scales= 'free_x' , 
		space = 'free_x' , 
	) + 
	ggplot2::scale_x_continuous(
		name= "Genomic Position (bp)", 
		breaks= seq( 0, max(toPlot$POS), 10^(floor(log10(max(toPlot$POS)))) ), 
		labels= format_genomic(), 
		expand= c(0,0)    # No expansion for plotting 
	) + 
	ggplot2::scale_y_continuous(
		name=vv, 
		expand= c(0,0)    # No expansion for plotting 
	) + 
	expand_limits( x=0, y=0 ) # Include the origin; 

p
dev.off()



