# [Sunhh@bioinfor01 fst]$ head -4 fstWC_CLAm_CLAs_w50ks5k.windowed.weir.fst
# CHROM   BIN_START       BIN_END N_VARIANTS      WEIGHTED_FST    MEAN_FST
# WM97pbV1_Chr01    1       50000   44      0.0272407       0.0193716
# WM97pbV1_Chr01    5001    55000   45      0.0275214       0.0204675
# WM97pbV1_Chr01    10001   60000   43      0.0267282       0.0192965

argvs <- commandArgs( trailingOnly=TRUE ) ;
if (is.na(argvs[1])) {
	message("\nRscript extract_top1_vcfFst_wind.R fstWC_CLAm_CLAs_w50ks5k.windowed.weir.fst outPrefix\n")
	q()
}
fn   <- argvs[1]   # fstWC_CLAm_CLAs_w50ks5k.windowed.weir.fst 
if (is.na(argvs[2])) {
	opref <- fn
} else {
	opref <- as.character( argvs[2] )
}
maxChrLen <- 40e6
if (is.na(argvs[3])) {
	; 
} else {
	maxChrLen <- as.numeric( argvs[3] )
}
chrlist <- NULL
if (is.na(argvs[4])) {
	; 
} else {
	chrlist <- read.table( as.character(argvs[4]), stringsAsFactors=FALSE )
	chrlist <- chrlist[,1]
}

fv <- read.table( file=fn, header=T, stringsAsFactors=F )
# Filter out organelles and NA sites. 
pdf( file=paste0(opref,'.hist.pdf') )
k_meanK <- fv$CHROM != 'WM97pbV1_Chr00' & fv$MEAN_FST >= 0
hist( fv$MEAN_FST[ k_meanK ] , breaks=seq(0,1,0.01), main=fn, xlab="Window's mean Fst" )
( f_mean_qt <- quantile( fv$MEAN_FST[ k_meanK ], probs=c(0,0.1,0.5,0.9,0.95,0.99,0.999,1) ) )
if (f_mean_qt[3] < 0.5) {
	l.lr <- 'topright'
} else {
	l.lr <- 'topleft'
}

legend(
	l.lr, 
	legend=c(
		paste0('Top   5%=', round(f_mean_qt[5], digits=4)),
		paste0('Top   1%=', round(f_mean_qt[6], digits=4)), 
		paste0('Top 0.1%=', round(f_mean_qt[7], digits=4)) 
	)
)
dev.off()

pp_fst <- function (cid = 'WM97pbV1_Chr01', fv=fv, ...) {
	fk <- fv$CHROM == cid
	mai <- par('mai')
	par('mai'=c(1.02, 0.82, 0.2, 0.42))
	plot( fv$BIN_START[fk], fv$MEAN_FST[fk], pch=20, ylim=c(0,1), xlab=cid, ylab="Fst", type='n', bty='l', xaxs='i', yaxs='i', ... )
	rect( 0, 0, max(fv$BIN_START[fk]), 1, density=NULL, col='grey', border='transparent' )
	points( fv$BIN_START[fk], fv$MEAN_FST[fk], pch=20, type='l' )
	segments( x0=0, y0=f_mean_qt[5], x1=max(fv$BIN_START[fk]), y1=f_mean_qt[5], col='purple'  )
	segments( x0=0, y0=f_mean_qt[6], x1=max(fv$BIN_START[fk]), y1=f_mean_qt[6], col='blue'  )
	segments( x0=0, y0=f_mean_qt[7], x1=max(fv$BIN_START[fk]), y1=f_mean_qt[7], col='red'  )
	par(mai=mai)
}

pdf( file=paste0(opref,'.chr.pdf'), height=3.5, width=21 )
if (is.null(chrlist)) {
	chrlist <- paste0('WM97pbV1_Chr', c(paste0('0',1:9),10,11))
}
# for (cid in paste0('WM97pbV1_Chr', c(paste0('0',1:9),10,11,'00'))) {
for (cid in chrlist) {
	pp_fst( cid=cid, fv=fv, xlim=c(0, maxChrLen) )
	legend( 'topright', legend=c(round(f_mean_qt[6], digits=4), round(f_mean_qt[7], digits=4)), col=c('blue', 'red'), pch=20 )
}
dev.off()

fv.top050 <- fv[ fv$MEAN_FST >= f_mean_qt[5], ]
fv.top010 <- fv[ fv$MEAN_FST >= f_mean_qt[6], ]
fv.top001 <- fv[ fv$MEAN_FST >= f_mean_qt[7], ]
write.table( fv.top050, file=paste0(opref,'.top0p050'), quote=FALSE, sep="\t", row.names=FALSE )
write.table( fv.top010, file=paste0(opref,'.top0p010'), quote=FALSE, sep="\t", row.names=FALSE )
write.table( fv.top001, file=paste0(opref,'.top0p001'), quote=FALSE, sep="\t", row.names=FALSE )

