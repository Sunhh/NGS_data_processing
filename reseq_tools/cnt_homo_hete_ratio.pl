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

Count genotyped (not missing) number, homozygous ratio and heteryzygous ratio per site. 
Out format : qw/chr pos GenoNum HomoRatio HeteRatio/

perl $0 in.snp_tbl > in.snp_tbl.cntHomHetR

-startColN       [$opts{'startColN'}]

Please note that the geno_col is $opts{'startColN'}

HH
!@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

unless ($opts{'noHeader'}) {
	my $l = <>; 
}
print join("\t",qw/chr pos GenoNum HomoRatio HeteRatio/)."\n"; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	# my $tot = $#ta- $geno_col +1; 
	my $tot = 0; 
	my $homN = 0; 
	my $hetN = 0; 
	for my $tb ( @ta[$geno_col .. $#ta] ) {
		($tb eq 'N' or $tb eq 'n') and next; 

		( $tb =~ m/^[ATGC*]$/i or $tb eq '*' or $tb =~ m/\+/ ) and do { $homN++; $tot++; next; }; 
		(&SNP_tbl::dna_d2b( &SNP_tbl::b2d($tb) )) > 1 and do { $hetN++; $tot++; next; }; 
		&tsmsg("[Wrn] Weired genotype [$tb] is treated as homozygous.\n"); 
		$homN++; $tot ++; 
	}
	my $ratHom = sprintf("%.4f", $homN/$tot) * 100; 
	my $ratHet = sprintf("%.4f", $hetN/$tot) * 100; 
	print join("\t", $ta[0], $ta[1], $tot, $ratHom, $ratHet)."\n"; 
}
