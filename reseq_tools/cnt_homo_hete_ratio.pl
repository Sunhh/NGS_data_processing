#!/usr/bin/perl 
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

my $geno_col = $opts{'startColN'}; 

sub usage {
	print STDERR <<HH; 

perl $0 in.snp_tbl > in.snp_tbl.cntHomHetR

-startColN       [$opts{'startColN'}]

Please note that the geno_col is $opts{'startColN'}
Out format : qw/ChromID Pos GenoN HomoRatio HeteRatio/

HH
	exit(1); 
}
!@ARGV and &usage(); 
$opts{'help'} and &usage(); 

my $l = <>; 
print join("\t",qw/ChromID Pos GenoN HomoRatio HeteRatio/)."\n"; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	# my $tot = $#ta- $geno_col +1; 
	my $tot = 0; 
	my $homN = 0; 
	my $hetN = 0; 
	for my $tb ( @ta[$geno_col .. $#ta] ) {
		if ( $tb =~ m/^[ATGC*]$|^[ATGC]\+[ATGC]+$/ ) {
			$homN ++; 
			$tot ++; 
		} elsif ( $tb eq 'N' ) {
			; 
		} elsif ( $tb =~ m/^[ATGC*][ATGC*]$/ ) {
			$hetN ++; 
			$tot ++; 
		} else {
			&tsmsg("Abnormal [$tb]\n"); 
		}
	}
	my $ratHom = sprintf("%.4f", $homN/$tot) * 100; 
	my $ratHet = sprintf("%.4f", $hetN/$tot) * 100; 
	print join("\t", $ta[0], $ta[1], $tot, $ratHom, $ratHet)."\n"; 
}
