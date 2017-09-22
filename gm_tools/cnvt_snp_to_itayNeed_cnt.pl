#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 

!@ARGV and die "perl $0 out_pref parent_LA483_UC204B_132offs.comm.snp\n"; 
my $pref = shift; 


# ==> parent_LA483_UC204B_132offs.comm.snp <==
# chr     pos     base    LA483   UC-204B Itay19  Itay2   Itay55  Itay75  RIL-1-G RIL-1-R RIL-11-G        RIL-11R RIL-12-G        RIL-12-R        RIL-13-G        RIL-13-
# SL2.40ch00      549274  T       G/G     T/T     T/T     T/T     T/T     G/G     T/T     T/T     T/T     T/G     ./.     T/T     ./.     T/T     T/T     T/T     T/T
# SL2.40ch00      551416  C       T/T     C/C     C/C     C/C     C/C     ./.     C/C     C/C     ./.     ./.     C/C     ./.     ./.     C/C     ./.     C/C     C/C
# SL2.40ch00      1146669 C       T/T     C/C     C/C     C/C     C/C     ./.     C/C     ./.     C/C     C/C     ./.     ./.     C/C     C/C     C/C     C/C     C/C


#my $hd = <>; 
#{
#	chomp($hd); 
#	my @ta=split(/\t/, $hd); 
#	my $b1 = "marker\tposition(bp)"; 
#	my $b2 = join("\t", @ta[5 .. $#ta]); 
#	$hd = join("\t", $b1, $b2)."\n"; 
#	print join("\t", @ta[0,1], qw/Total A_N B_N AF%/)."\n"; 
	print join("\t", qw/chr pos/, qw/Total A_N B_N AF%/)."\n"; 
#}

my %have; 

while (<>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	my $b1 = "$ta[0]_$ta[1]"; 
	my $fn = "${pref}$ta[0].csv"; 
	if ($. == 1 and $ta[0] eq "chr") {
		next; 
	}
#	unless (defined $have{$fn}) {
#		$have{$fn} = 1; 
#		&fileSunhh::write2file($fn, $hd, '>' ); 
#	}
	my $g1 = $ta[3]; 
	my $g2 = $ta[4]; 
	$g1 =~ s!^([ATGC])/\1$!$1! or die "g1=$g1\n"; 
	$g2 =~ s!^([ATGC])/\1$!$1! or die "g2=$g2\n"; 
	my %cnt; 
	for my $tc (@ta[ 5 .. $#ta ]) {
		if ( $tc eq "./." ) {
			$tc = "-"; 
		} elsif ( $tc eq "$g1/$g1" ) {
			$tc = "a"; 
		} elsif ( $tc eq "$g2/$g2" ) {
			$tc = "b"; 
		} elsif ( $tc eq "$g1/$g2" or $tc eq "$g2/$g1" ) {
			$tc = "h"; 
		} else {
			# die "tc=|$tc| g1=$g1 g2=$g2\n"; 
			$tc = '-'; 
		}
		$cnt{$tc} ++; 
	}
	delete($cnt{"-"}); 
	my $af = 'NA'; 
	for my $a (qw/a b h/) {
		$cnt{$a} //= 0; 
	}
	my $ttl = 2*$cnt{'a'} + 2*$cnt{'b'} + 2*$cnt{'h'}; 
	my $ca = $cnt{'a'} * 2 + $cnt{'h'}; 
	my $cb = $cnt{'b'} * 2 + $cnt{'h'}; 
	$ttl > 0 and $af = sprintf("%.2f", $ca / $ttl * 100); 
	print join("\t", @ta[0,1], $ttl, $ca, $cb, $af)."\n"; 
	# &fileSunhh::write2file($fn, join("\t", $b1, @ta[5 .. $#ta])."\n", '>>'); 
}

