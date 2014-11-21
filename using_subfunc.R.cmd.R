####### clean Q20 PE

source("/home/Sunhh/tools/clean_reads/using_subfunc.R")
pattern_lis <- read.table("pattern_lis.03", header=T, stringsAsFactors=F)
# pattern_lis <- read.table("pattern_lis.02", header=T, stringsAsFactors=F)
# pattern_lis <- read.table("pattern_lis.01", header=T, stringsAsFactors=F)

for ( i in 1:nrow(pattern_lis) ) {
inFq1 <- paste0(pattern_lis$Prefix[i], "_R1.ndupB", sep="")
inFq2 <- paste0(pattern_lis$Prefix[i], "_R2.ndupB", sep="")
oFq1 <- paste0(pattern_lis$Prefix[i], "_R1", sep="")
oFq2 <- paste0(pattern_lis$Prefix[i], "_R2", sep="")
adp1 <- pattern_lis$R1pattern[i]
adp2 <- pattern_lis$R2pattern[i]
clean.pe.fq.file( inFqName1=inFq1, outFqName1=oFq1, adaptor1=adp1, inFqName2=inFq2, outFqName2=oFq2, adaptor2=adp2, RdPerYield=10e6 )
}

####### clean MP
source("/home/Sunhh/tools/clean_reads/using_subfunc.R")  # Penguin server.
# source("/home/Sunhh/tools/github/data_proc/using_subfunc.R")  # WWZ server. 

pattern_lis <- read.table("MP_junc_pattern.01", header=T, stringsAsFactors=F)
# pattern_lis <- read.table("MP_junc_pattern.02", header=T, stringsAsFactors=F)
# pattern_lis <- read.table("MP_junc_pattern.03", header=T, stringsAsFactors=F)
for ( i in 1:nrow(pattern_lis) ) {
gc()
inFq1 <- paste0(pattern_lis$Prefix[i], "_R1.paired", sep="")
inFq2 <- paste0(pattern_lis$Prefix[i], "_R2.paired", sep="")
oFq1 <- paste0(pattern_lis$Prefix[i], sep="")
junc_seq <- pattern_lis$JuncPattern[i]
	myseq <- junc_seq
	myseq_comp <- chartr("ATGC", "TACG", myseq)
	substring(myseq[1], 1:nchar(myseq[1]), 1:nchar(myseq[1]))
	x <- strsplit(myseq_comp, "")
	x <- lapply(x, rev)
	myseq_revcomp <- sapply(x, paste, collapse="")
junc_seq <- c(junc_seq, myseq_revcomp)
clean.mp.fq.file ( inFqName1=inFq1, outFqName1=oFq1, inFqName2=inFq2, junction.seq=junc_seq, RdPerYield=20e6, min.length=40, align.opts=list() )
}

