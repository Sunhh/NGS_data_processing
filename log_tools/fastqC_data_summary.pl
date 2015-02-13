#!/usr/bin/perl -w 
# Used to extract reads number information from fastqC results. 
use strict; 

!@ARGV and die "perl $0 *_fastqc/fastqc_data.txt\n"; 

my %dataInfo; 
for my $file (@ARGV) {
	open F,'<',"$file" or die "file=$file\n"; 
	while (<F>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		if ($ta[0] eq 'Filename') {
			$dataInfo{$file}{fn} = $ta[1]; 
		} elsif ($ta[0] eq 'File type') {
			$dataInfo{$file}{fileType} = $ta[1]; 
		} elsif ($ta[0] eq 'Encoding') {
			$dataInfo{$file}{Encoding} = $ta[1]; 
		} elsif ($ta[0] eq 'Total Sequences') {
			$dataInfo{$file}{seqN} = $ta[1]; 
		} elsif ($ta[0] eq 'Filtered Sequences' or $ta[0] eq 'Sequences flagged as poor quality') {
			$dataInfo{$file}{filtN} = $ta[1]; 
		} elsif ($ta[0] eq 'Sequence length') {
			$dataInfo{$file}{seqLen} = $ta[1]; 
		} elsif ($ta[0] eq '%GC') {
			$dataInfo{$file}{GCperc} = $ta[1]; 
		}
	}
	close F; 
	print STDOUT join("\t", $file, $dataInfo{$file}{fn}, $dataInfo{$file}{Encoding}, $dataInfo{$file}{seqN}, $dataInfo{$file}{seqLen}, $dataInfo{$file}{GCperc}, $dataInfo{$file}{filtN}, $dataInfo{$file}{fileType})."\n"; 
}

# ##FastQC        0.10.1
# >>Basic Statistics      pass
# #Measure        Value
# Filename        P1_200d_R2.fastq.gz
# File type       Conventional base calls
# Encoding        Sanger / Illumina 1.9
# Total Sequences 46833520
# Filtered Sequences      0
# Sequence length 101
# %GC     49
#

