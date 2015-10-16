#!/usr/bin/perl 
use strict; 
use warnings; 
use SNP_tbl; 
my $st_obj = SNP_tbl->new(); 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"startColN:i", # 2 
	"noHeader!", 
); 
$opts{'startColN'} //= 2; 

my $geno_col = $opts{'startColN'}; 

my $help_txt = <<HH; 

perl $0 in.snp_tbl > in.snp_tbl.cntHomHetR

-startColN       [$opts{'startColN'}]

Please note that the geno_col is $opts{'startColN'}
Out format : qw/ChromID Pos GenoN HomoRatio HeteRatio/

HH
!@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

unless ($opts{'noHeader'}) {
	my $l = <>; 
}
print join("\t",qw/chr pos GenoN HomoRatio HeteRatio/)."\n"; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	# my $tot = $#ta- $geno_col +1; 
	my $tot = 0; 
	my $homN = 0; 
	my $hetN = 0; 
	for my $tb ( @ta[$geno_col .. $#ta] ) {
		($tb eq 'N' or $tb eq 'n') and next; 

		( $ta[$i] =~ m/^[ATGC*]$/i or $ta[$i] eq '*' or $ta[$i] =~ m/\+/ ) and do { $homN++; $tot++; next; }; 
		(&SNP_tbl::dna_d2b($ta[$i])) > 1 and do { $hetN++; $tot++; next; }; 
		&tsmsg("[Wrn] Weired genotype [$tb] is treated as homozygous.\n"); 

		$tot ++; 
		
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
