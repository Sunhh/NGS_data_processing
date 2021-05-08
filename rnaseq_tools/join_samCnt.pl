#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"toJoinFile:s@", 
	"toJoinHead:s@", 
	"toJoinList:s", 
	"skipCntLine:i", 
	"useEleIDList:s", 
	"help!", 
); 

my %gg; 
&set_glob(); 

my @outTbl; 
if ( defined $gg{'useEleID'} and @{$gg{'useEleID'}} > 0 ) {
	push(@outTbl, ['EleID']); 
	for (@{$gg{'useEleID'}}) {
		push(@outTbl, [$_]); 
	}
}

for (my $i=0; $i<@{$gg{'bam2id'}}; $i++) {
	my ($cntFn, $sampleID) = @{$gg{'bam2id'}[$i]}; 
	my $cntFh = &openFH($cntFn, '<'); 
	my $out_idx0 = 0; 
	if ( defined $gg{'useEleID'} and @{$gg{'useEleID'}} > 0 ) {
		push(@{$outTbl[0]}, $sampleID); 
		my %id2v; 
		while (<$cntFh>) {
			$. <= $gg{'skipCntLine'} and next; 
			chomp; 
			my @ta = split(/\t/, $_); 
			$id2v{$ta[0]} //= $ta[1]; 
		}
		for my $a1 (@outTbl) {
			my $v = (defined $id2v{$a1->[0]}) ? $id2v{$a1->[0]} : 'NA' ; 
			push(@$a1, $v); 
		}
	} else {
		if ($i == 0) {
			$outTbl[0] = ['EleID', $sampleID]; 
			while (<$cntFh>) {
				$. <= $gg{'skipCntLine'} and next; 
				# $out_idx0 ++; 
				chomp; 
				my @ta = split(/\t/, $_); 
				push(@outTbl, [@ta[0,1]]); 
			}
		} else {
			push(@{$outTbl[0]}, $sampleID); 
			while (<$cntFh>) {
				$. <= $gg{'skipCntLine'} and next; 
				$out_idx0 ++; 
				chomp; 
				my @ta = split(/\t/, $_); 
				$outTbl[$out_idx0][0] eq $ta[0] or &stopErr("[Err] Mismatch between index ID and following ID. [$outTbl[$out_idx0][0] VS. $ta[0]] in file $cntFn\n"); 
				push(@{$outTbl[$out_idx0]}, $ta[1]); 
			}
		}
	}
	close($cntFh); 
}
for (@outTbl) {
	print STDOUT join("\t", @$_)."\n"; 
}

&tsmsg("[Rec] $0 done.\n"); 

########################################################
# Sub-routines 
########################################################

sub set_glob {
	$gg{'skipCntLine'} = 0; 
	$gg{'help_txt'} = <<"H1"; 
################################################################################
# perl $0 -toJoinList samCntFiles_to_paste
#
# -help
#
# -toJoinList            [filename] Format : path_to_bamCnt \\t bamCnt_name
#
# -toJoinFile            [filename] Paired with -toJoinHead ; path_to_bamCnt ; multiple times; 
# -toJoinHead            [filename] Paired with -toJoinFile ; bamCnt_name    ; multiple times; 
#
# -skipCntLine           [number]  default '$gg{'skipCntLine'}'. No. of lines skipped from combination; 
#
# -useEleIDList          [filename] A file name recording Elements ID I want to extract from cnt files. 
#                                   This will replace the order and content in real bamCnt files. 
################################################################################
H1
	$opts{'help'} and &LogInforSunhh::usage($gg{'help_txt'}); 
	my @bam2id; 
	if ( defined $opts{'toJoinFile'} ) {
		if ( defined $opts{'toJoinHead'} ) {
			@{$opts{'toJoinFile'}} == @{$opts{'toJoinHead'}} or &stopErr ("[Err] Please pair -toJoinFile and -toJoinHead files. "); 
		} else {
			@{$opts{'toJoinHead'}} = @{$opts{'toJoinFile'}}; 
		}
		for (my $i=0; $i<@{$opts{'toJoinFile'}}; $i++) {
			push(@bam2id, [ $opts{'toJoinFile'}[$i], $opts{'toJoinHead'}[$i] ]); 
		}
	}
	if ( defined $opts{'toJoinList'} ) {
		my $fh = &openFH($opts{'toJoinList'}, '<'); 
		while (<$fh>) {
			m!^\s*(#|$)! and next; 
			chomp($_); 
			my @ta = &splitL("\t", $_); 
			$ta[1] //= $ta[0]; 
			push(@bam2id, [$ta[0], $ta[1]]); 
		}
		close($fh); 
	}
	@bam2id == 0 and &stopErr("[Err] Please assign ' -toJoinList ' or ' -toJoinFile ' \n"); 
	$gg{'bam2id'} = \@bam2id; 

	defined $opts{'skipCntLine'} and $gg{'skipCntLine'} = $opts{'skipCntLine'}; 

	if (defined $opts{'useEleIDList'}) {
		my $fh = &openFH($opts{'useEleIDList'}); 
		while (<$fh>) {
			chomp; 
			my @ta = split(/\t/, $_); 
			push(@{$gg{'useEleID'}}, $ta[0]); 
		}
		close($fh); 
	}
	return; 
}# set_glob() 


