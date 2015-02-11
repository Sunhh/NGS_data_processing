source("./using_subfunc.R")
pattern_lis <- read.table("MP_junc_pattern", header=T, stringsAsFactors=F)
for ( i in 1:nrow(pattern_lis) ) {
gc()
inFq1 <- paste0(pattern_lis$Prefix[i], "_R1.ndupB", sep="")
inFq2 <- paste0(pattern_lis$Prefix[i], "_R2.ndupB", sep="")
oFq1  <- paste0(pattern_lis$Prefix[i], sep="")
junc_seq <- pattern_lis$JuncPattern[i]
	myseq <- junc_seq
	myseq_comp <- chartr("ATGC", "TACG", myseq)
	substring(myseq[1], 1:nchar(myseq[1]), 1:nchar(myseq[1]))
	x <- strsplit(myseq_comp, "")
	x <- lapply(x, rev)
	myseq_revcomp <- sapply(x, paste, collapse="")
junc_seq <- c(junc_seq, myseq_revcomp)

clean.mp.fq.file ( inFqName1=inFq1, outFqName1=oFq1, inFqName2=inFq2, junction.seq=junc_seq, RdPerYield=5e6, min.length=40, max.chunks=10)
}

