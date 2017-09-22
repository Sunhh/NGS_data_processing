# ==> set01_13_to_11.fst.perSiteChrPos <==
# chr     pos     Ho      Hs      Ht      Dst     Htp     Dstp    Fst     Fstp    Fis     Dest
# chr10   467     0       0       0       0       0       0       NA      NA      NA      0

argvs <- commandArgs( trailingOnly=TRUE ) ;
fn   <- argvs[1]   # in.fst.perSiteChrPos
topR <- as.numeric( argvs[2] ) # 0.01 
if (topR > 1) {
	print (topR)
	quit()
}

aa <- read.table( file=fn, header=T, stringsAsFactors=F )
# Filter out organelles and NA sites. 
aa.kk <- aa$chr != "plast" & aa$chr != "mito" & !is.na(aa$Fst) & aa$Fst >= 0
aa <- aa[ aa.kk, ]
aa.kk <- NULL 

# Find the threshold of topR 
aa.qt <- quantile( aa$Fst, probs=c(0,0.5,0.95,1-topR,1) )
aa.thres <- aa.qt[4]
cat("threshold for", topR, "is ", aa.thres, "\n")

# Get the selected sites. 
aa.slct <- aa[ aa$Fst >= aa.thres, ]
write.table( aa.slct, file=paste0( fn, ".top", topR, sep=""), append=F, row.names=F, col.names=T, quote=F, sep="\t" )


