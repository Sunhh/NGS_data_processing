#!/usr/bin/perl
# http://gtrnadb.ucsc.edu/faq.html
# The tRNA candidates that meet the following two criteria are considered as potential tRNA-derived SINEs and are eliminated from the final prediction results displayed in the GtRNAdb: 
### (1) a final bit score of less than 55 bits (50 bits for tRNA-SeC-TCA)
### (2) are identified by only one of the two pre-filter scanners (trnascan 1.4 and EufindtRNA)
### As such, the tRNA genes listed within the GtRNAdb will differ from stand-alone tRNAscan-SE analyses due to this ad-hoc post-filtering step. The full set of unfiltered tRNAs + likely tRNA-like SINES can be provided upon request.

# Accoding to tRNAscan-SE v2.0 result (using Infernal) and the Ath data from GtRNAdb reference on 201806, 
#   the above two filters have not been applied to the dataset. 
# And the direct output of tRNAscan-SE v2.0 is quite similar to the reference : 
#   1) There are 642 tRNAs predicted by v2.0 program, while 641 in the reference; 
#   2) There are 19 tRNAs from program which are not exactly same to the reference, and in these 19 tRNAs, 
#   2.1) 15 tRNAs are mostly overlapped with the reference, having final bit scores ranging from 20.4 to 71.0; 
#   2.2) Only four tRNAs don't overlap the reference, having low final bit scores (Max=34.1). 
#   3) There are 18 tRNAs from the reference which are not exactly same to the program, and in them, 
#   3.1) The three non-overlapping tRNAs have very low final bit score (Max=26.7); 
# So for the purpose of predicting tRNAs, I want to accept the default parameter of tRNAscan-SE v2.0, and do not do any furhter filtering after that. 



use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 

!@ARGV and die "perl $0 arab.nucl.fasta.trnase.result arab.nucl.fasta.trnase.resultLegacy > arab.nucl.fasta.trnase.result.filt\n"; 
my $ifile_1 = shift; 
my $ifile_2 = shift; 
my %support2; 
if (defined $ifile_2 and $ifile_2 ne '') {
	my $ifh_2 = &openFH($ifile_2, '<'); 
	while (<$ifh_2>) {
		chomp; 
		my @ta = &splitL("\t", $_); 
		if ( !(defined $ta[1]) or $ta[1] !~ m!^\d+\s*$! ) {
			next; 
		}
		# (2) are identified by only one of the two pre-filter scanners (trnascan 1.4 and EufindtRNA)
		$ta[11] eq 'Bo' or next; 
		$ta[0] =~ s!^\s+|\s+$!!g; 
		$ta[2] =~ s!^\s+|\s+$!!g; 
		$ta[3] =~ s!^\s+|\s+$!!g; 
		my $tk = join("\t", $ta[0], $ta[2], $ta[3]); 
		$support2{$tk} = 1; 
	}
	close($ifh_2); 
}
# Input file 'arab.nucl.fasta.trnase.resultLegacy' comes from command of tRNAscan-SE v2.0 with legacy method: 
#  tRNAscan-SE \ 
#  -E -y -L --detail --progress --forceow -H \
#  -o arab.nucl.fasta.trnase.resultLegacy \
#  arab.nucl.fasta
# [Sunhh@bioinfor01 aa]$ head -5 arab.nucl.fasta.trnase.resultL
# Sequence                tRNA            Bounds          tRNA    Anti    Intron Bounds   Cove    HMM     2'Str   Hit
# Name            tRNA #  Begin           End             Type    Codon   Begin   End     Score   Score   Score   Origin  Note
# --------        ------  -----           ------          ----    -----   -----   ----    ------  -----   -----   ------  ------
# Chr1            1       306384          306456          Val     TAC     0       0       79.63   51.80   27.83   Bo
# Chr1            2       515494          515566          Phe     GAA     0       0       77.43   54.21   23.22   Bo
# Chr1            7       1159023         1159093         Gly     GCC     0       0       57.43   45.84   11.59   Eu
# Chr1            8       1324923         1325007         Met     CAT     1324961 1324971 60.94   40.83   20.11   Bo
# Chr1            15      2751329         2751401         Arg     CCG     0       0       37.42   36.58   0.84    Eu      pseudo
#

# Input file 'arab.nucl.fasta.trnase.result' comes from command of tRNAscan-SE v2.0 with Infernal method: 
#  tRNAscan-SE \
#   -E -y --detail --progress --forceow -H \
#   -o arab.nucl.fasta.trnase.result \
#   -m# -f# -l# -b# -a# \
#   arab.nucl.fasta
#
# [Sunhh@bioinfor01 aa]$ head arab.nuclsta.out
# Sequence                tRNA            Bounds          tRNA    Anti    Intron Bounds   Inf     HMM     2'Str   Hit     Isotype Isotype
# Name            tRNA #  Begin           End             Type    Codon   Begin   End     Score   Score   Score   Origin  CM      Score   Note
# --------        ------  -----           ------          ----    -----   -----   ----    ------  -----   -----   ------  ------- ------- ------
# Chr1            1       306384          306456          Val     TAC     0       0       74.4    44.00   30.40   Inf     Val     103.3
# Chr1            2       515494          515566          Phe     GAA     0       0       78.3    58.40   19.90   Inf     Phe     113.2
# Chr1            3       552640          552711          His     GTG     0       0       63.0    33.70   29.30   Inf     His     86.8
# Chr1            4       604402          604474          Lys     CTT     0       0       87.7    67.60   20.10   Inf     Lys     106.7

my $ifh_1 = &openFH($ifile_1, '<'); 
while (<$ifh_1>) {
	chomp; 
	my @ta = &splitL("\t", $_); 
	if ( !(defined $ta[1]) or $ta[1] !~ m!^\d+\s*$! ) {
		print STDOUT "$_\n"; 
		next; 
	}
	$ta[0] =~ s!^\s+|\s+$!!g; 
	$ta[2] =~ s!^\s+|\s+$!!g; 
	$ta[3] =~ s!^\s+|\s+$!!g; 
	$ta[4] =~ s!^\s+|\s+$!!g; 
	$ta[5] =~ s!^\s+|\s+$!!g; 
	$ta[8] =~ s!^\s+|\s+$!!g; 
	# (1) a final bit score of less than 55 bits (50 bits for tRNA-SeC-TCA)
	if ($ta[4] eq 'SeC') {
		$ta[5] eq 'TCA' or &stopErr("[Err] Bad line : $_\nta[4]=$ta[4] ta[5]=$ta[5]\n"); 
		$ta[8] >= 50 or next; 
	} else {
		$ta[8] >= 55 or next; 
	}
	# (2) supported by two legacy methods; 
	if ( defined $ifile_2 and $ifile_2 ne '' ) {
		my $tk = join("\t", $ta[0], $ta[2], $ta[3]); 
		defined $support2{$tk} or next; 
	}
	print STDOUT "$_\n"; 
}



