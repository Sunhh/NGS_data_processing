#!/usr/bin/perl -w
# Used to check duplicated reads in sweetpotato reads.
# Edit for WM_D1_test 
# Edit 2012-10-09
# Edit 2013-08-12 Add a parameter to cut sequences as key. 
# Edit 2013-09-30 Add -subseqS to give start position to cut read sequences. (1-index based)
# Editing 2014-03-03 Add -maxDupNum to give a maximum number of duplicates output. Add -rcDup to calculate duplicates with reverse complementary ones. 
use strict;
use Getopt::Long; 


my $ptag = 'out'; 
my $help = 0; 
my $is_ogz = 0; 
my $no_out = 0; 
my $subseqP = 0; 
my $subseqS = 1; 
my $rcDup = 0; 
my $maxDupN = 1; 
GetOptions(
'opre=s' => \$ptag, 
'help' => \$help, 
'ogz' => \$is_ogz,
'noout' => \$no_out, 
'subseq=i'  => \$subseqP, 
'subseqS=i' => \$subseqS, 
'rcDup' => \$rcDup, 
'maxDupNum=i' => \$maxDupN
); 

sub usage {
	my $version = 'v1.0'; 
	my $last_time = '2014-03-03'; 
	print STDOUT <<HELP;
#******* Instruction of this program *********#
Version   : $version
Last time : $last_time
perl $0 -opre out_prefix input1.fastq input2.fastq

-help     
-ogz          [Boolean] output gzipped files. Default OFF. 
-noout        [Boolean] No output files. Only calculate duplicates. Default OFF. 
-subseqS      [INT] Start position from which to pick subsequences. default 1. 
-subseq       [INT] Length of subsequences. Default 0 for whole length. 
-rcDup        [Boolean] Remove paired-duplicates considering reverse complement fragment. Default OFF. 
-maxDupNum    [INT] Maximum number of copies kept in .ndupB files. Default 1. 
#******* Instruction of this program *********#

HELP
	exit(1); 
}

( !@ARGV or $help ) and &usage(); 
( -f $ARGV[0] and -f $ARGV[1] ) or do { &timeLog("\n[Err] Input files [@ARGV] not found.\n\n"); &usage(); }; 

## # Do not change 1-based start position to 0-based position for substr() function. 
# $subseqS--; 
$subseqS < 1 and $subseqS = 1; 
$subseqP < 0 and $subseqP = 0; # to ignore extract sequence from the tail. 

my ($rf1, $rf2) = @ARGV; 

if ($rf1 =~ m/\.gz$/) {
	open (R1,'-|', "gzip -cd $rf1") or die; 
}else{
	open (R1, '<', "$rf1") or die; 
}
if ($rf2 =~ m/\.gz$/) {
	open (R2,'-|', "gzip -cd $rf2") or die; 
}else{
	open (R2, '<', "$rf2") or die; 
}

&timeLog( "[Rec] Reading $ptag files: $rf1 $rf2\n" ); 
&timeLog( "[Rec] Parameters: -subseqS=$subseqS -subseq=$subseqP -rcDup=$rcDup -maxDupNum=$maxDupN -noout=$no_out -opre=$ptag\n" ); 

unless ($no_out) {
	if ($is_ogz) {
		open (O1,'|-', "gzip -c >${ptag}_R1.ndupB.gz") or die;
		open (O2,'|-', "gzip -c >${ptag}_R2.ndupB.gz") or die; 
		open (D1,'|-', "gzip -c >${ptag}_R1.dropB.gz") or die; 
		open (D2,'|-', "gzip -c >${ptag}_R2.dropB.gz") or die;
		&timeLog( "[Rec] Create output files: ${ptag}_R1.ndupB.gz ${ptag}_R2.ndupB.gz ${ptag}_R1.dropB.gz ${ptag}_R2.dropB.gz\n" ); 
	}else{
		open (O1,'>', "${ptag}_R1.ndupB") or die;
		open (O2,'>', "${ptag}_R2.ndupB") or die;
		open (D1,'>', "${ptag}_R1.dropB") or die; 
		open (D2,'>', "${ptag}_R2.dropB") or die; 
		&timeLog( "[Rec] Create output files: ${ptag}_R1.ndupB ${ptag}_R2.ndupB ${ptag}_R1.dropB ${ptag}_R2.dropB\n" ); 
	}
}

my (%have1, %have2); 
my %have_both; 
my $total = 0; 
my $keepN = 0; 
my $dropN = 0; 
my $uniq = 0; 
my $mult = 0; 
while (!eof(R1)) {
	$. % 10000000 == 0 and &timeLog( "[Msg] $. line in $ptag files. Kept $keepN read pairs.\n" ); 
	my ($k1, $s1, $q1) = &readArecord(\*R1);
	my ($k2, $s2, $q2) = &readArecord(\*R2);
	$total ++; 
	my $os1 = $s1; 
	my $os2 = $s2; 
	unless ($subseqP == 0) {
		$s1 = substr($s1, $subseqS-1, $subseqP); 
		$s2 = substr($s2, $subseqS-1, $subseqP); 
	}
	my $s_both = "$s1\t$s2"; 
	$have_both{$s_both} ++; 
	my $dupN; 
	if ( $rcDup ) {
		my $s_both_rc = "$s2\t$s1"; 
		$have_both{$s_both_rc} ++; 
		$dupN = ($have_both{$s_both}+$have_both{$s_both_rc})/2; 
	} else {
		$dupN = $have_both{$s_both}; 
	}

	if ( $dupN == 2 ) {
		$uniq --; 
		$mult ++; 
	}
	if ($dupN > $maxDupN) {
		$dropN ++; 
		unless ( $no_out ) {
			print D1 "$k1\n$os1\n+\n$q1\n"; 
			print D2 "$k2\n$os2\n+\n$q2\n"; 
		}
	} else {
		unless ( $no_out ) {
			print O1 "$k1\n$os1\n+\n$q1\n"; 
			print O2 "$k2\n$os2\n+\n$q2\n"; 
		}
		$keepN ++; 
		$dupN == 1 and $uniq ++; 
	}
}

unless ($no_out) {
	close D1; 
	close D2; 
	close O1;
	close O2;
}
close R1;
close R2;
my $uniq_R = &rate( $uniq, $total ); 
my $mult_R = &rate( $mult, $total ); 
my $drop_R = &rate( $dropN, $total ); 
my $keep_R = &rate( $keepN, $total ); 

&timeLog( "[Rec] Finish $ptag files.\n" ); 
&timeLog( "[Rec] There are $total read pairs in total [$ptag].\n" ); 
&timeLog( "[Rec] There are $uniq ($uniq_R\%) unique patterns, $mult ($mult_R\%) multiple patterns, and $keepN ($keep_R\%) reads kept in both (subseq=$subseqP subseqS=$subseqS).\n" ); 

sub readArecord {
		my ($fh) = @_;
		my $tk = <$fh>;
		my $ts = <$fh>;
		<$fh>;
		my $tq = <$fh>;
		$tk =~ s/\s+$//; 
		$ts =~ s/\s+$//;
		$tq =~ s/\s+$//;
		return ($tk, $ts, $tq);
}

sub rate {
	my ($n1, $n2) = @_; 
	if ($n2 == 0) {
		warn "Denominator is zero.\n"; 
		return 'NA'; 
	}else{
		my $r = int($n1/$n2*10000+0.5)/100; 
		return $r; 
	}
}

sub timeLog {
	my $tt = scalar(localtime()); 
	print STDERR join('', "[$tt] ", @_); 
}


