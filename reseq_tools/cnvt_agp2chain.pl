#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 

my $help_info = <<HH; 
######################################################################
# perl $0 ctg.fas.key_len chr.fas.key_len ctg2chr.agp > ctg2chr.agp.chain
######################################################################
#  By : sunhonghe_1984@163.com
#  This is used to convert AGP file to UCSC chain format for picard-liftovervcf function. 
#  Usage example : 
#   perl $0 wm97pb_v2ID.scf.fa.kl WM97pbV1.chr.fa.kl WM97pbV1.scf2chr.agp WM97pbV1.scf2chr.agp.chain
#   java -jar /home/Sunhh/src/align/picard/v2.16.0/picard/build/libs/picard.jar CreateSequenceDictionary  R=WM97pbV1.chr.fa  O=WM97pbV1.chr.dict
#   java -jar /home/Sunhh/src/align/picard/v2.16.0/picard/build/libs/picard.jar LiftoverVcf \
#     LIFTOVER_MIN_MATCH=0 \
#     WARN_ON_MISSING_CONTIG=ture \
#     WRITE_ORIGINAL_POSITION=ture \
#     I=refPBv2_pass.vcf \
#     O=refPBv2_pass_scf2Chr.vcf \
#     CHAIN=WM97pbV1.scf2chr.agp.chain \
#     REJECT=refPBv2_pass_scf2Chr_rejected.vcf \
#     R=WM97pbV1.chr.fa
#### 
#  Try sepRun perl for speed and memory concerns. 
######################################################################
HH

@ARGV >= 2 or &LogInforSunhh::usage($help_info); 

!@ARGV and -t and die "perl $0 ctg_keyLen chr_keyLen input.agp > output.chain\n"; 

# https://genome.ucsc.edu/goldenPath/help/chain.html
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_vcf_LiftoverVcf.php

# CHR : target / reference; 
# CTG : query 

#      chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id 
#      chain 4900 chrY 58368225 + 25985406 25985566 chr5 151006098 - 43549808 43549970 2
#      16      0       2
#      60      4       0
#      10      0       4
#      70

my %chrLen; 
my %ctgLen; 
my $ctgF = shift; 
%ctgLen = map { $_->[0] => $_->[1] } &fileSunhh::load_tabFile( $ctgF ); 
my $chrF = shift; 
%chrLen = map { $_->[0] => $_->[1] } &fileSunhh::load_tabFile( $chrF ); 

my @out; 
my $chain_n = 0; 
while (<>) {
	m!^\s*#! and next; 
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[4] eq 'N' and next; 
	$ta[4] eq 'W' or die "ta[4]=$ta[4]: $_\n"; 
	$chain_n ++; 
	defined $chrLen{$ta[0]} or die "No length for chr [$ta[0]]\n"; 
	defined $ctgLen{$ta[5]} or die "No length for ctg [$ta[5]]\n"; 
	if ( $ta[8] =~ m!^(\+|\?)$! ) {
		#                        chain   score tName   tSize            tStrand tStart    tEnd    qName   qSize            qStrand qStart    qEnd    id 
		print STDOUT join("\t", 'chain', 100,  $ta[5], $ctgLen{$ta[5]}, '+',    $ta[6]-1, $ta[7], 
			                               $ta[0], $chrLen{$ta[0]}, '+',    $ta[1]-1, $ta[2], 
			                $chain_n
		)."\n"; 
	} elsif ( $ta[8] =~ m!^(\-)$! ) {
		#                        chain   score tName   tSize            tStrand tStart                  tEnd    qName   qSize            qStrand qStart                  qEnd                      id 
		print STDOUT join("\t", 'chain', 100,  $ta[5], $ctgLen{$ta[5]}, '+',    $ta[6]-1, $ta[7], 
			                               $ta[0], $chrLen{$ta[0]}, '-',    $chrLen{$ta[0]}-$ta[2], $chrLen{$ta[0]}-$ta[1]+1, 
			                $chain_n
		)."\n"; 
	}
	print STDOUT ($ta[2]-$ta[1]+1)."\n"; 
	print STDOUT "\n"; 
}


