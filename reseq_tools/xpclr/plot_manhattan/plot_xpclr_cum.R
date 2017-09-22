#!/home/Sunhh/bin/Rscript

.tsmsg <- function(...) {
	message("[", date(), "]: ", ...)
}# End .tsmsg

argvs <- commandArgs( trailingOnly=TRUE ) ;

fn_jn     <- as.character( argvs[1] )
fn_cumChr <- as.character( argvs[2] )
fn_outpdf <- paste0(fn_jn, '.mht_plot.pdf')
fn_jn_02  <- NULL
if ( !is.na(argvs[3]) ) { fn_outpdf <- as.character( argvs[3] ) }
if ( !is.na(argvs[4]) ) { fn_jn_02  <- as.character( argvs[4] ) }
if ( length(argvs) < 2 ) {
	.tsmsg("Rscript plot_xpclr_cum.R xpclr_w10ks10k chrLen_cum [out_pdf xpclr_w10ks10k_rev]\n"); 
	q()
}

main_txt    <- c('', '')
add_mainTxt <- FALSE
add_mainTxt <- TRUE
if ( add_mainTxt ) {
	main_txt[1] <- fn_jn
	if (!is.null(fn_jn_02)) { main_txt[2] <- fn_jn_02 }
}

plot_jnChr <- function ( jn=jn , chrlen=chrlen , cv='Avg', cols=NULL, jn_02=NULL, main_txt=c('', ''), ... ) {
	if ( is.null(cols) ) { cols <- rainbow(length(chrlen$chrID)) }
	k1 <- !is.na( jn[[cv]] )
	y_max <- max( jn[[cv]][k1] )
	y_min <- 0
	x_max <- max( jn$WindE[k1] )
	# x_min <- min( jn$WindS )
	x_min <- 0
	k2 <- NULL 
	if ( !is.null(jn_02) ) {
		k2 <- !is.na( jn_02[[cv]] )
		y_max <- max( jn[[cv]][k1] , jn_02[[cv]][k2] )
		x_max <- max( jn$WindE[k1] , jn_02$WindE[k2] )
		# x_min <- min( jn$WindS , jn_02$WindS )
		par( mfrow=c(2,1) )
	}
	plot( c(x_min, x_max), c(y_min, y_max), type='n', axes=FALSE , xlab=NA , ylab='XP-CLR Score', main=main_txt[1], ... )
	axis( side=1, at = c(chrlen$chrCumS, chrlen$chrCumE) , labels = FALSE , tick = TRUE )
	axis( side=1, at = rowMeans( chrlen[,c(3,4)]) , labels = chrlen$chrID , tick = FALSE )
	axis( side=2 )
	for (i in 1:length(chrlen$chrID)) {
		kk <- jn$WindS[k1] >= chrlen$chrCumS[i] & jn$WindE[k1] <= chrlen$chrCumE[i]
		points( jn$WindS[k1][kk] , jn[[cv]][k1][kk] , pch=1, cex=0.5, col=cols[i], ... )
		# points( jn$WindS[kk] , jn[[cv]][kk] , pch='.',           col=cols[i] )
	}
	if ( !is.null(jn_02) ) {
		plot( c(x_min, x_max), c(y_min, y_max), type='n', axes=FALSE , xlab=NA , ylab='XP-CLR Score', main=main_txt[2], ... )
		axis( side=1, at = c(chrlen$chrCumS, chrlen$chrCumE) , labels = FALSE , tick = TRUE )
		axis( side=1, at = rowMeans( chrlen[,c(3,4)]) , labels = chrlen$chrID , tick = FALSE )
		axis( side=2 )
		for (i in 1:length(chrlen$chrID)) {
			kk <- jn_02$WindS[k2] >= chrlen$chrCumS[i] & jn_02$WindE[k2] <= chrlen$chrCumE[i]
			points( jn_02$WindS[k2][kk] , jn_02[[cv]][k2][kk] , pch=2, cex=0.5, col=cols[i], ... )
			# points( jn_02$WindS[kk] , jn_02[[cv]][kk] , pch='.',           col=cols[i] )
		}
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
data_jn_02  <- NULL 
if ( !is.null(fn_jn_02) ) {
	data_jn_02 <- read.table( fn_jn_02, header=T, stringsAsFactors=F )
}

pdf( file=fn_outpdf , height=7, width=21 )
plot_jnChr( jn= data_jn, chrlen= data_cumChr, cv='Avg', jn_02=data_jn_02, main_txt=main_txt )
dev.off()


