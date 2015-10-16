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
	"acceptHete!", 
); 
$opts{'startColN'} //= 2; 
use SNP_tbl; 
my $st_obj = SNP_tbl->new(); 

my $geno_col = $opts{'startColN'}; 
my $help_txt = <<HH; 

perl $0 in.snp > in.snp.diff_mat

-help
-startColN      [$opts{'startColN'}]
-acceptHete     [Bool] Accepting heterozygous sites, and compare them like this : 
                  'AG' vs. 'AG' accounts 0   for diff. 
                  'AG' vs. 'AA' accounts 0.5 for diff. 
                  'AG' vs. 'TT' accounts 1   for diff. 
                  'AG' vs. 'TA' accounts 0.5 for diff. 


HH
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
	if ( $opts{'acceptHete'} ) {
		for (my $i=$geno_col; $i<@hh; $i++) {
			my $h1 = &sep_allele( $ta[$i] ); 
			my @a1 = keys %$h1; 
			( my $n1 = scalar(@a1) ) > 0 or next; 
			$n1 <= 2 or &stopErr("Why here! [n1=$n1, \@a1=@a1]\n"); 
			#my $b1 = &fmt_base( $ta[$i] ); 
			#$b1 eq 'N' and next; 
			for (my $j=$i+1; $j<@hh; $j++) {
				my $h2 = &sep_allele( $ta[$j] ); 
				my @a2 = keys %$h2; 
				( my $n2 = scalar(@a2) ) > 0 or next; 
				$n2 <= 2 or &stopErr("Why here! [n2=$n2, \@a2=@a2]\n"); 
				#my $b2 = &fmt_base( $ta[$j] ); 
				#$b2 eq 'N' and next; 
				my $kk = "${i}-${j}"; 
				$cnt{$kk}{'total'} ++; 
				
				if ($n1 == 1) {
					if ( $n2 == 1 ) {
						$a1[0] eq $a2[0] or $cnt{$kk}{'diffN'} ++; 
					} else {
						$cnt{$kk}{'diffN'} += ( ( defined $h2->{$a1[0]} ) ? 0.5 : 1 ); 
					}
				} else {
					for my $a1_b ( @a1 ) {
						defined $h2->{$a1_b} or $cnt{$kk}{'diffN'} += 0.5; 
					}
				}
				#$b1 eq $b2 or $cnt{$kk}{'diffN'} ++; 
			}
		}
	} else {
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
	$a =~ m/\*|\+/ and return $a; 
	$a = $st_obj->SingleChar($a, 'onlyATGC'=>1 ); 
	return $a; 
}

sub sep_allele {
	# Only two alleles allowed. 
	my %h; 
	if ($_[0] =~ m/\*|\+/) {
		$h{$_[0]} = 1; 
	} else {
		my $t1 = $st_obj->SingleChar( $_[0], 'maxAlleleN'=>2 ); 
		$t1 eq 'N' and return \%h; 
		for ( &SNP_tbl::dna_d2b($t1) ) {
			$h{$_} = 1; 
		}
	}
	return (\%h); 
}

