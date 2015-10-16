#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use SNP_tbl; 
my $st_obj = SNP_tbl->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"startColN:i", # 2 
	"noHeader!", 
); 
$opts{'startColN'} //= 2; 

my $type_startColN = $opts{'startColN'}; 

my $help_txt = <<HH; 

Count N (missing), heterozygous, and homozygous sites per individual. 
The output format is : qw/IndvID N_Num Typed_Num Het_Num Hom_Num Het_Ratio Hom_Ratio N_Ratio/

perl $0 in.snp > in.snp.cntNHH

-help               
-startColN          [$opts{'startColN'}]
-noHeader           [Bool]

Genotype column start from colN=$opts{'startColN'}
Do not parse the first line. 

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my @ha; 
unless ($opts{'noHeader'}) {
	my $head = <>; 
	chomp($head); 
	@ha=split(/\t/, $head); 
}
my @cnt_N; 
my @cnt_hete; 
my @cnt_homo; 

while (<>) { 
	($.-1) % 1e6 == 1 and &tsmsg("[Msg] $. lines.\n"); 
	chomp; 
	my @ta = split(/\t/, $_); 
	unless ($opts{'noHeader'} and @ha == 0) {
		for (my $i=0; $i<@ta; $i++) {
			$ha[$i] = "Col_$i"; 
		}
	}
	for (my $i=$type_startColN; $i<@ta; $i++) { 
		($ta[$i] eq 'N' or $ta[$i] eq 'n') and do { $cnt_N[$i]++; next; }; 
		( $ta[$i] =~ m/^[ATGC*]$/i or $ta[$i] eq '*' or $ta[$i] =~ m/\+/ ) and do { $cnt_homo[$i]++; next; }; # The '*' and sites with \+ are treated as homozygous. I don't like genotype like 'A*', so I want to remove it before calculation. 
		(&SNP_tbl::dna_d2b( &SNP_tbl::b2d($ta[$i]) )) > 1 and do { $cnt_hete[$i]++; next; }; # Heterozygous. 
		&tsmsg("[Wrn] Weired genotype [$ta[$i]] is treated as homozygous.\n"); 
		$cnt_homo[$i] ++; 
		# $ta[$i]=~m/^[ATGC*]$|^[ATGC]\+[ATGC]+$/ and $cnt[$i]++; 
	} 
} 
print STDOUT join("\t", qw/IndvID N_Num Typed_Num Het_Num Hom_Num Het_Ratio Hom_Ratio N_Ratio/)."\n"; 
for (my $i=0; $i<@ha; $i++) { 
	$cnt_N[$i] //= 0; 
	$cnt_hete[$i] //= 0; 
	$cnt_homo[$i] //= 0; 
	my $tot = $cnt_N[$i]+$cnt_hete[$i]+$cnt_homo[$i]; 
	my $tot_typed = $cnt_hete[$i]+$cnt_homo[$i]; 
	my $rat_hete = ($tot_typed > 0) ? sprintf("%02.02f", $cnt_hete[$i]/$tot_typed*100) : -1 ; 
	my $rat_homo = ($tot_typed > 0) ? sprintf("%02.02f", $cnt_homo[$i]/$tot_typed*100) : -1 ; 
	my $rat_miss = ($tot       > 0) ? sprintf("%02.02f", $cnt_N[$i]   /$tot_typed*100) : -1 ; 
	print "$ha[$i]\t$cnt_N[$i]\t$tot_typed\t$cnt_hete[$i]\t$cnt_homo[$i]\t$rat_hete\t$rat_homo\t$rat_miss\n";
}


