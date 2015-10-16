#!/usr/bin/perl -w
# Rules: 
#  1. Skip 'N'; 
#  2. Skip '+/*'; 
#  3. Only consider homo genotype. 
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"startColN:i", # 2 
); 
$opts{'startColN'} //= 2; 
use SNP_tbl; 
my $st_obj = SNP_tbl->new(); 

my $geno_col = $opts{'startColN'}; 
sub usage {
	print STDERR <<HH;

perl $0 in.snp > in.snp.diff_mat

-help
-startColN      [$opts{'startColN'}]

\$geno_col is $geno_col


HH
	exit(1); 
}

!@ARGV and &usage(); 
$opts{'help'} and &usage(); 

my @hh; 
my %cnt; 
while (<>) {
	chomp; 
	@hh = split(/\t/, $_); 
	last; 
}

while (<>) {
	$. % 1000 == 1 and &tsmsg("[Msg] $. line.\n"); 
	chomp; 
	my @ta = split(/\t/, $_); 
	for (my $i=$geno_col; $i<@hh; $i++) {
		my $b1 = &fmt_base( $ta[$i] ); 
		$b1 eq 'N' and next; 
		for (my $j=$i+1; $j<@hh; $j++) {
			my $b2 = &fmt_base( $ta[$j] ); 
			$b2 eq 'N' and next; 
			my $kk = "${i}-${j}"; 
			$cnt{$kk}{'total'} ++; 
			$b1 eq $b2 or $cnt{$kk}{'diffN'} ++; 
		}
	}
}

print STDOUT join("\t", qw/Indv1 Indv2 NumTotal NumDiff/)."\n"; 
for my $kk ( sort { my ($a1,$a2) = ( $a =~ m!^(\d+)\-(\d+)$! ); my ($b1,$b2) = ( $b =~ m!^(\d+)\-(\d+)$! ); $a1 <=> $b1 || $a2 <=> $b2;  } keys %cnt ) {
	my ($k1, $k2) = ( $kk =~ m!^(\d+)\-(\d+)$! ); 
	$cnt{$kk}{'total'} //= 0; 
	$cnt{$kk}{'diffN'} //= 0; 
	print STDOUT join("\t", $hh[$k1], $hh[$k2], $cnt{$kk}{'total'}, $cnt{$kk}{'diffN'})."\n"; 
}

sub fmt_base {
	my $a = shift; 
	$a = uc($a); 
	$a = $st_obj->SingleChar($a, 'onlyATGC'=>1 ); 
	return $a; 
}


