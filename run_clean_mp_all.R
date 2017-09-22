# using_subfunc.R

source("./using_subfunc.R")

# Remove junction sequences from Mate-paired 
MP_pattern_lis <- read.table("MP_junc_pattern", header=T, stringsAsFactors=F)
for ( i in 1:nrow(MP_pattern_lis) ) {
gc()
inFq1 <- paste0(MP_pattern_lis$Prefix[i], "_R1.ndupB", sep="")
inFq2 <- paste0(MP_pattern_lis$Prefix[i], "_R2.ndupB", sep="")
oFq1  <- paste0(MP_pattern_lis$Prefix[i], sep="")
junc_seq <- MP_pattern_lis$JuncPattern[i]
	myseq <- junc_seq
	myseq_comp <- chartr("ATGC", "TACG", myseq)
	substring(myseq[1], 1:nchar(myseq[1]), 1:nchar(myseq[1]))
	x <- strsplit(myseq_comp, "")
	x <- lapply(x, rev)
	myseq_revcomp <- sapply(x, paste, collapse="")
junc_seq <- c(junc_seq, myseq_revcomp)

clean.mp.fq.file ( inFqName1=inFq1, outFqName1=oFq1, inFqName2=inFq2, junction.seq=junc_seq, RdPerYield=10e6, min.length=40, max.chunks=20)
}


# Remove low quality and adaters from separated Mate-paired reads. 
PE_pattern_lis <- read.table("PE_pattern_lis", header=T)
for ( i in 1:nrow(PE_pattern_lis) ) {
inFq1 <- paste0(PE_pattern_lis$Prefix[i], ".p1a", sep="")
inFq2 <- paste0(PE_pattern_lis$Prefix[i], ".p2a", sep="")
oFq1 <- paste0(PE_pattern_lis$Prefix[i], "_a_R1", sep="")
oFq2 <- paste0(PE_pattern_lis$Prefix[i], "_a_R2", sep="")
adp1 <- PE_pattern_lis$R1pattern[i]
adp2 <- PE_pattern_lis$R2pattern[i]
clean.pe.fq.file( inFqName1=inFq1, outFqName1=oFq1, adaptor1=adp1, inFqName2=inFq2, outFqName2=oFq2, adaptor2=adp2, RdPerYield=10e6 )

inFq1 <- paste0(PE_pattern_lis$Prefix[i], ".p1b", sep="")
inFq2 <- paste0(PE_pattern_lis$Prefix[i], ".p2b", sep="")
oFq1 <- paste0(PE_pattern_lis$Prefix[i], "_b_R1", sep="")
oFq2 <- paste0(PE_pattern_lis$Prefix[i], "_b_R2", sep="")
adp1 <- PE_pattern_lis$R1pattern[i]
adp2 <- PE_pattern_lis$R2pattern[i]
clean.pe.fq.file( inFqName1=inFq1, outFqName1=oFq1, adaptor1=adp1, inFqName2=inFq2, outFqName2=oFq2, adaptor2=adp2, RdPerYield=10e6 )

inFq1 <- paste0(PE_pattern_lis$Prefix[i], ".s1", sep="")
oFq1 <- paste0(PE_pattern_lis$Prefix[i], "_s_R1", sep="")
adp1 <- PE_pattern_lis$R1pattern[i]
clean.pe.fq.file( inFqName1=inFq1, outFqName1=oFq1, adaptor1=adp1, RdPerYield=10e6 )

inFq1 <- paste0(PE_pattern_lis$Prefix[i], ".s2", sep="")
oFq1 <- paste0(PE_pattern_lis$Prefix[i], "_s_R2", sep="")
adp1 <- PE_pattern_lis$R1pattern[i]
clean.pe.fq.file( inFqName1=inFq1, outFqName1=oFq1, adaptor1=adp1, RdPerYield=10e6 )
}

