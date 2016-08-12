#!/home/Sunhh/bin/Rscript

.tsmsg <- function(...) {
	message("[", date(), "]: ", ...)
}# End .tsmsg

argvs <- commandArgs( trailingOnly=TRUE ) ;

fn_jn     <- as.character( argvs[1] )
fn_cumChr <- as.character( argvs[2] )
if ( length(argvs) == 0 ) {
	.tsmsg("Rscript plot_xpclr_cum.R xpclr_w10ks10k chrLen_cum\n"); 
	q()
}

plot_jnChr <- function ( jn=jn , chrlen=chrlen , cv='Avg', cols=NULL, ... ) {
	if ( is.null(cols) ) { cols <- rainbow(length(chrlen$chrID)) }
	plot( c(min(jn$WindS), max(jn$WindE)), c(0, max(jn[[cv]])), type='n', axes=FALSE , xlab=NA , ylab='XP-CLR Score', ... )
	axis( side=1, at = c(chrlen$chrCumS, chrlen$chrCumE) , labels = FALSE , tick = TRUE )
	axis( side=1, at = rowMeans( chrlen[,c(3,4)]) , labels = chrlen$chrID , tick = FALSE )
	axis( side=2 )
	for (i in 1:length(chrlen$chrID)) {
		kk <- jn$WindS >= chrlen$chrCumS[i] & jn$WindE <= chrlen$chrCumE[i]
		points( jn$WindS[kk] , jn[[cv]][kk] , pch=20, cex=0.3, col=cols[i], ... )
		# points( jn$WindS[kk] , jn[[cv]][kk] , pch='.',           col=cols[i] )
	}
}
plot_sepChr_slct <- function ( sep_all=sep_all , sep_slct=sep_slct, chrID='chr1', ... ) {
	k_all  <- sep_all[,1]  == chrID & !is.na( sep_all[,6] )
	k_slct <- sep_slct[,1] == chrID & !is.na( sep_slct[,6] )
	plot(   c(-1, max(sep_all[k_all,2])), c(-1, max(sep_all[k_all,6])), type='n', xlab=chrID, ylab='XP-CLR Score', ... )
	points( sep_all[k_all,2], sep_all[k_all,6], pch=20, cex=1, col='blue', ... )
	segments( x0=sep_slct[k_slct, 2], y0= -1, x1=sep_slct[k_slct, 3], y1= -1, col='red', lwd=5, ... )
}

data_jn     <- read.table( fn_jn,     header=T, stringsAsFactors=F )
data_cumChr <- read.table( fn_cumChr, header=T, stringsAsFactors=F )

pdf( file=paste0(fn_jn, '.jnChr.pdf') , height=7, width=21 )
plot_jnChr( jn= data_jn, chrlen= data_cumChr, cv='Avg' )
dev.off()


