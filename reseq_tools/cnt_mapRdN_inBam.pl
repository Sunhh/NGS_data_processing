#!/usr/bin/perl
use strict; 
use warnings; 
use SeqAlnSunhh; 
use LogInforSunhh; 

my $htxt = <<HHH; 
################################################################################
perl $0 in.bam > in.bam.map_stat

Result (in.bam.map_stat) example: 
Filename                  All          ProperPairedRdN   PairedRdN    SingleRdN     UnmappedRdN    OtherRdN
F582R1_bySM_fix.bam       57827298     55968230          1542390      108379        208296         3
F582R2HG_bySM_fix.bam     56482040     54124958          1905874      141387        309792         29
SQ026_bySM_fix.bam        285572114    276544304         4574622      397791        4055112        285
################################################################################
HHH

!@ARGV and die "$htxt\n"; 

my %flag_unmap = %{ &SeqAlnSunhh::mk_flag( 'keep' => '2=1' ) }; 
my %flag_r1    = %{ &SeqAlnSunhh::mk_flag( 'keep' => '6=1,7=0' ) };
my %flag_r2    = %{ &SeqAlnSunhh::mk_flag( 'keep' => '6=0,7=1' ) }; 
my %flag_goodPair = %{ &SeqAlnSunhh::mk_flag( 'keep' => '0=1,1=1,2=0,3=0,8=0,11=0' ) }; 
my %flag_errPair  = %{ &SeqAlnSunhh::mk_flag( 'keep' => '0=1,1=0,2=0,3=0,8=0,11=0' ) }; 
my %flag_1map     = %{ &SeqAlnSunhh::mk_flag( 'keep' => '0=1,1=0,2=0,3=1,8=0,11=0' ) }; 

# Setup values; 
my %scores; 
for (keys %flag_goodPair) { $scores{$_} = 1; }
for (keys %flag_errPair)  { $scores{$_} = 2; }
for (keys %flag_1map)     { $scores{$_} = 3; }
for (keys %flag_unmap)    { $scores{$_} = 4; }

my (%rd1,%rd2); 


my $fn = shift; 

open F,'-|',"samtools view $fn" or die ; 
while (<F>) {
	$. % 100e3 == 1 and &tsmsg("[Msg] $. line in [$fn]\n"); 
	chomp; 
	m!^(\S+)\t(\d+)\t! or die "$_\n"; 
	my ($rdid, $ff) = ($1, $2); 
	# defined $scores{$ff} or next; 
	defined $scores{$ff} or $scores{$ff} = 5; 
	if (defined $flag_r1{$ff}) {
		if (defined $rd1{$rdid}) {
			if ($rd1{$rdid} > $scores{$ff}) {
				$rd1{$rdid} = $scores{$ff}; 
			}
		} else {
			$rd1{$rdid} = $scores{$ff}; 
		}
	} elsif (defined $flag_r2{$ff}) {
		if (defined $rd2{$rdid}) {
			if ($rd2{$rdid} > $scores{$ff}) {
				$rd2{$rdid} = $scores{$ff}; 
			}
		} else {
			$rd2{$rdid} = $scores{$ff}; 
		}
	} else {
		die "r12: $_\n"; 
	}
}
close F; 

warn "a1\n"; 

my %cnt; 
for my $rdid (keys %rd1) {
	( defined $rd1{$rdid} and defined $rd2{$rdid} ) or die "rdid=$rdid\n"; 
	if ( $rd1{$rdid} eq $rd2{$rdid} ) {
		$cnt{ $rd1{$rdid} } += 2; 
	} else {
		$cnt{ $rd1{$rdid} } += 1; 
		$cnt{ $rd2{$rdid} } += 1; 
	}
}
warn "a2\n"; 
for (qw/1 2 3 4 5/) {
	$cnt{$_} //= 0; 
	$cnt{'total'} += $cnt{$_}; 
}

print join("\t", qw/Filename  All ProperPairedRdN  PairedRdN  SingleRdN  UnmappedRdN OtherRdN/)."\n"; 
print join("\t", $fn, @cnt{qw/total 1 2 3 4 5/})."\n"; 


