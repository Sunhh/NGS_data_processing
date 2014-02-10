fn <- 'NSP_20kb_GCCAAT_L00M.SP_gt25k.bam.ins'

fnames <- read.table("list", header=F, stringsAsFactors=F)

pdf(file="aa.pdf")
for (ii in 1:length(fnames$V1)) {
fn <- fnames$V1[ii]
l.ins = scan(fn)
( q.ins = quantile(l.ins, probs=c(0,0.01,0.05,0.10,0.25,0.5,0.75,0.90,0.95,0.99,1) ))
minIns <- 0
maxIns <- max( q.ins[10], 30000 )

stepIns <- 10 
use.ins <- l.ins >= minIns & l.ins <= maxIns
( mean.ins = mean(l.ins[use.ins]) )
( sd.ins = sd(l.ins[use.ins]) )
( median.ins = median(l.ins[use.ins]) )
# kk <- l.ins >= 0 & l.ins <= 20000
kk = use.ins
# pdf(file=paste0(fn,"_histo.pdf"))
( xmin <- minIns )
# xmax <- maxIns
( xmax <- 600 )
# ( xmax <- 12e3 )
# ( xmax <- 25e3 )
t.l.ins <- l.ins
t.l.ins[t.l.ins < xmin] <- xmin
t.l.ins[t.l.ins > xmax] <- xmax
hist( t.l.ins, breaks=seq(xmin,max(t.l.ins+stepIns), by=stepIns), main=fn, xlab="Insert length (bp)", xlim=c(xmin, xmax) )
# hist( t.l.ins, breaks=seq(xmin,max(t.l.ins+stepIns), by=stepIns), main=fn, xlab="Insert length (bp)", xlim=c(xmin, xmax), ylim=c(0,15) )
ppp <- ifelse( median.ins <= mean(c(xmin, xmax)) , 'topright', 'topleft')
legend(ppp,
legend=c(
  paste0("Median=", median.ins),
  paste0("Mean=", mean.ins),
  paste0("SD=", sd.ins),
  paste0("Max=", max(l.ins)), 
	paste0("Num=", length(l.ins))
)
)
# dev.off()
}
dev.off()


