
date
#cp SlidingWindowTrimmer.java org/usadellab/trimmomatic/trim/SlidingWindowTrimmer.java
#javac org/usadellab/trimmomatic/trim/SlidingWindowTrimmer.java
java org.usadellab.trimmomatic.Trimmomatic SE -phred33 CG_Levi_15kb_R1.ndupB CG_Levi_15kb_R1.ndupB.clean SLIDINGWINDOW:4:20
date


# test code

cp FastqRecord.java org/usadellab/trimmomatic/fastq/FastqRecord.java
cp SlidingWindowTrimmer.java org/usadellab/trimmomatic/trim/SlidingWindowTrimmer.java
javac org/usadellab/trimmomatic/trim/SlidingWindowTrimmer.java
javac org/usadellab/trimmomatic/fastq/FastqRecord.java
java org.usadellab.trimmomatic.Trimmomatic SE -phred33 test.fastq test.clean.fastq SLIDINGWINDOW:4:20
