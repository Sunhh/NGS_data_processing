#!/home/Sunhh/bin/Rscript

argvs <- commandArgs( trailingOnly=TRUE ) ;
# in.bam.ins , out_hist.pdf , plot_title , minIns , maxIns , stepIns ; 
if ( is.na(argvs[1]) ) {
	message("Rscript in.bam.ins out_hist.pdf plot_title minIns maxIns stepIns skipSmall skipBig")
	q()
}

fn <- as.character( argvs[1] ) ; 
outPdf <- paste0( fn, "_hist.pdf", sep="" ); 
if ( !is.na(argvs[2]) ) {
	outPdf <- as.character( argvs[2] ) 
}
title <- fn
if ( !is.na(argvs[3]) ) {
	title <- as.character( argvs[3] )
}

minIns <- -1
maxIns <- -1
if ( !is.na(argvs[4]) ) {
	minIns <- as.numeric(argvs[4])
}
if ( !is.na(argvs[5]) ) {
	maxIns <- as.numeric(argvs[5])
}
stepIns <- 1
if ( !is.na(argvs[6]) ) {
	stepIns <- as.numeric(argvs[6])
}
skipSmall <- FALSE
if ( !is.na(argvs[7]) ) {
	skipSmall <- as.logical( argvs[7] )
	if ( is.na(skipSmall)) {
		skipSmall <- FALSE 
		message("Bad input for skipSmall")
	}
}
skipBig <- FALSE
if ( !is.na(argvs[8]) ) {
	skipBig <- as.logical( argvs[8] )
	if ( is.na(skipBig) ) {
		skipBig <- FALSE
		message("Bad input for skipBig")
	}
}

l.ins <- scan( fn )
( q.ins = quantile(l.ins, probs=c(0,0.01,0.05,0.10,0.25,0.5,0.75,0.90,0.95,0.99,1) ))
if (minIns == -1) {
	minIns <- min(l.ins)
}
if (maxIns == -1) {
	maxIns <- max(l.ins)
}

cnt.ins <- l.ins >= q.ins[3] & l.ins <= q.ins[9]
( mean.ins = round(mean(l.ins[cnt.ins]), digits=2) )
( median.ins = median(l.ins[cnt.ins]) )
( sd.ins = round(sd(l.ins[cnt.ins]), digits=2) )

# Plot 
( xmin <- minIns )
( xmax <- maxIns )
t.l.ins <- l.ins
if ( skipSmall ) {
	t.l.ins <- l.ins[ l.ins >= xmin ]
} else {
	t.l.ins[t.l.ins < xmin] <- xmin
}
if ( skipBig ) {
	t.l.ins <- l.ins[ l.ins <= xmax ]
} else {
	t.l.ins[t.l.ins > xmax] <- xmax
}

pdf( file=outPdf )
hist( t.l.ins, breaks=seq(xmin,max(t.l.ins+stepIns), by=stepIns), main=title, xlab="Insert length (bp)", xlim=c(xmin, xmax) )
ppp <- ifelse( median.ins <= mean(c(xmin, xmax)) , 'topright', 'topleft')
legend(ppp,
legend=c(
  paste0("Median=", median.ins),
  paste0("Mean=", mean.ins),
  paste0("SD=", sd.ins),
#  paste0("Max=", max(l.ins)),
  paste0("Pairs=", length(l.ins))
)
)
# abline( v= mean.ins, col='red' )
dev.off()
message( paste0("Title=", title, " Mean=", mean.ins, "Median=", median.ins, "SD=", sd.ins, sep="") )

