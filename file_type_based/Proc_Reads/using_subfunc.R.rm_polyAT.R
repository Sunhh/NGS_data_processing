####### clean PolyAT SE 

source("./using_subfunc.R")
pref_lis <- read.table("pref_list", header=F, stringsAsFactors=F)

for ( i in 1:nrow(pref_lis) ) {
	inFq1 <- paste0( pref_lis$V1[i], '_R1.fq' , sep='' )
	oFq1  <- paste0( pref_lis$V1[i], '_trimAT.fq', sep='' )
	adp1  <- 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
		myseq <- adp1
		myseq_comp <- chartr('ATGC', 'TACG', myseq) 
		substr( myseq[1], 1:nchar(myseq[1]), 1:nchar(myseq[1]) ) 
		x <- strsplit( myseq_comp, '' )
		x <- lapply(x, rev) 
		myseq_revcomp <- sapply(x, paste, collapse='') 
	adp1 <- c( adp1, myseq_revcomp )
	clean.pe.fq.file( inFqName1= inFq1, outFqName1= oFq1, adaptor1= adp1, RdPerYield= 10e6 , qual.opts=list( min.qual=0 ) )
}

