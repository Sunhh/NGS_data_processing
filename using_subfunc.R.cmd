
# Record: 
#  I found in P1_200a_R1.fastq.gz and P1_300_ATCACG_R1.fastq.gz, there are R1 reads end with 'ATCTCGTATGCCGTCTTCTGCTTG' without barcode or 'AGATCGGAAGAGCACACGTCTGAACTiCCAGTCAC' sequences. So for these cases, I might need to run two times of function-clean.pe.fq.file with different "adaptor1" option. 
# And ~10bp head of 'ATCTCGTATGCCGTCTTCTGCTTG' is similar to some chloroplast sequence or rRNA sequence. Like this 'ATGTGTGTCGGTTGCTTGCA'; 
# I can also find 'ATGTGTGTCGGTTGCTTGCA' in P1_300 library, maybe it is common in P1 libraries. So I plan to trim twice and see how many reads changed. 
# Not sure if it is common to other species. 

# Need to see the read infor: 
# [Sunhh@Penguin SeqData]$ gzip -cd P1_500b_R1.fastq.gz | grep ATCACGATCTCGTATGCCGTCTTCTGCTTG
# CAATCTTAACATCAAGGTTGAGATCGGAAGAGCACACGCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAACAGAATTAGAATCAT



# For PE adaptor cleaning. 

barcod <- 'AAATTT'
patternR1 <- paste0( "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", barcod, "ATCTCGTATGCCGTCTTCTGCTTG",sep="")
patternR2 <- "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"

barcode_lis <- read.table("barcode_lis", header=T, stringsAsFactors=F)
pattern_lis <- NULL
for ( i in 1:nrow(barcode_lis) ) {
	t.prefix <- barcode_lis$Prefix[i]
	t.R1pattern <- paste0( "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", barcode_lis$Barcode[i], "ATCTCGTATGCCGTCTTCTGCTTG",sep="" )
	t.R2pattern <- "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
	pattern_lis <- rbind( t.prefix, barcode_lis$Barcode[i], t.R1pattern, t.R2pattern )
}
colnames(pattern_lis) <- c("Prefix", "Barcode", "R1pattern", "R2pattern")

source("/home/Sunhh/tools/github/data_proc/using_subfunc.R")
pattern_lis <- read.table("pattern_lis", header=T, stringsAsFactors=F)
for ( i in 1:nrow(pattern_lis) ) {
inFq1 <- paste0(pattern_lis$Prefix[i], "_R1.ndupB", sep="")
inFq2 <- paste0(pattern_lis$Prefix[i], "_R2.ndupB", sep="")
oFq1 <- paste0(pattern_lis$Prefix[i], "_R1", sep="")
oFq2 <- paste0(pattern_lis$Prefix[i], "_R2", sep="")
adp1 <- pattern_lis$R1pattern[i]
adp2 <- pattern_lis$R2pattern[i]
clean.pe.fq.file( inFqName1=inFq1, outFqName1=oFq1, adaptor1=adp1, inFqName2=inFq2, outFqName2=oFq2, adaptor2=adp2, RdPerYield=40e6 )

system(paste( "perl /home/Sunhh/tools/fq_rdNum.pl ",
        oFq1, ".paired.hq ",
        oFq2, ".paired.hq ",
        oFq1, ".single.hq ",
        oFq2, ".single.hq ",
        oFq1, ".paired ",
        oFq2, ".paired ",
        oFq1, ".single ",
        oFq2, ".single ",
        " >> test_tbl ; ",

        "rm -f ",
        oFq1, ".paired.hq ",
        oFq2, ".paired.hq ",
        oFq1, ".single.hq ",
        oFq2, ".single.hq ",

        sep="")
)
}

clean.mp.fq.file( inFqName1=inFq1, outFqName1=oFq1, inFqName2=inFq2, junction.seq=junc_seq, RdPerYield=40e5 )


# [Sunhh@Penguin Q20]$ more MP_junc_pattern
Prefix  JuncPattern
P1_cre5k_CACTCA_time1   CGTATAACTTCGTATAATGTATGCTATACGAAGTTATACA
P1_cre5k_CACTCA_time2   CGTATAACTTCGTATAATGTATGCTATACGAAGTTATACA
P1_ec5k_ATGAGC_time1    CGTATAACTTCGTATAATGTATGCTATACGAAGTTATACA
P1_ec5k_ATGAGC_time2    CGTATAACTTCGTATAATGTATGCTATACGAAGTTATACA
P1_MP2k GGTCGATAACTTCGTATAATGTATGCTATACGAAGTTATACA
P1_MP5k GGTCGATAACTTCGTATAATGTATGCTATACGAAGTTATACA

source("/home/Sunhh/tools/clean_reads/using_subfunc.R")  # Penguin server.
pattern_lis <- read.table("MP_junc_pattern", header=T, stringsAsFactors=F)
for ( i in 1:nrow(pattern_lis) ) {
gc()
inFq1 <- paste0(pattern_lis$Prefix[i], "_R1.paired", sep="")
inFq2 <- paste0(pattern_lis$Prefix[i], "_R2.paired", sep="")
oFq1 <- paste0(pattern_lis$Prefix[i], sep="")
junc_seq <- pattern_lis$JuncPattern[i]
clean.mp.fq.file ( inFqName1=inFq1, outFqName1=oFq1, inFqName2=inFq2, junction.seq=junc_seq, RdPerYield=40e6 )
}


