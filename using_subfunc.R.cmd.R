### Change trimmomatic directory for different server !!! 

#### For PE clean from .ndupB 
source( '/home/Sunhh/tools/clean_reads/using_subfunc.R' ); 
pattern_lis <- read.table( 'pattern_lis.01', header=T, stringsAsFactors=F ); 
# pattern_lis <- read.table( 'pattern_lis.02', header=T, stringsAsFactors=F ); 
# pattern_lis example : 
# Prefix  Barcode R1pattern                                                               R2pattern
# P3_200  ATCACG  AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG        AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
# P3_500a CCGTCC  AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCATCTCGTATGCCGTCTTCTGCTTG        AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
# P3_500b GAGTGG  AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATCTCGTATGCCGTCTTCTGCTTG        AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
for ( i in 1:nrow(pattern_lis) ) {
	inFq1 <- paste0( pattern_lis$Prefix[i], "_R1.ndupB", sep="" ); 
	inFq2 <- paste0( pattern_lis$Prefix[i], "_R2.ndupB", sep="" ); 
	oFq1  <- paste0( pattern_lis$Prefix[i], "_R1", sep="" ); 
	oFq2  <- paste0( pattern_lis$Prefix[i], "_R2", sep="" ); 
	clean.pe.fq.file( inFqName1=inFq1, outFqName1=oFq1, adaptor1=adp1, inFqName2=inFq2, outFqName2=oFq2, adaptor2=adp2, RdPerYield=10e6 ); 
}


#### For MP clean from .ndupB ####
source("./using_subfunc.R")
pattern_lis <- read.table("MP_junc_pattern", header=T, stringsAsFactors=F)
# MP_junc_pattern example : 
#   If processing PE together, should add columns [ 'R1pattern' , 'R2pattern' ] for PE. 
# Prefix          JuncPattern             R1pattern                                                               R2pattern
# PI482246_5kbA   CTGTCTCTTATACACATCT     AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG        GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
# PI482246_20kb   CTGTCTCTTATACACATCT     AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG      GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
# PI482246_200A   N                       AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCATCTCGTATGCCGTCTTCTGCTTG        GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
# PI482246_500B   N                       AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAATCTCGTATGCCGTCTTCTGCTTG        GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

for ( i in 1:nrow(pattern_lis) ) {
	gc(); 
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

	clean.mp.fq.file ( inFqName1=inFq1, outFqName1=oFq1, inFqName2=inFq2, junction.seq=junc_seq, RdPerYield=5e6, min.length=40, max.chunks=10, align.opts=list() )


	gc(); 
	adp1  <- pattern_lis$R1pattern[i]; 
	adp2  <- pattern_lis$R2pattern[i]; 
	for ( j in c("a", "b") ) {
		inFq1_j <- paste0( oFq1, ".p1", j, sep="" ); 
		inFq2_j <- paste0( oFq1, ".p2", j, sep="" ); 
		oFq1_j  <- paste0( oFq1, "_", j, "_Q20_R1", sep="" ); 
		oFq2_j  <- paste0( oFq1, "_", j, "_Q20_R2", sep="" ); 
		
		clean.pe.fq.file( inFqName1=inFq1_j, outFqName1=oFq1_j, adaptor1=adp1, inFqName2=inFq2_j, outFqName2=oFq2_j, adaptor2=adp2, RdPerYield=10e6, max.chunks=10 )
	}
}

