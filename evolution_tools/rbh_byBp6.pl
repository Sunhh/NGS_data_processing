#!/usr/bin/perl
use strict; 
use warnings; 
use SeqAlnSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
use fileSunhh; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"in_bp6:s@", 

	"min_similarity:f", 
	"max_lenDiffR:f", 
	"log_lines:i", # 0
	"gene1_list:s", 
	"gene2_list:s", 
	"combine_gene12!", 
	"paired_need!", 
#	"min_cscore:f", 
); 

$opts{'min_similarity'} //= 0; 
$opts{'max_lenDiffR'} //= 0; 
#$opts{'min_cscore'} //= 0; 
$opts{'log_lines'} //= 0; 

my $help_txt = <<HH; 

perl $0 -in_bp6 A2B_B2A.bp6 > A2B_B2A.bp6.rbh.pairs
OR
perl $0 -in_bp6 A2B.bp6 -in_bp6 B2A.bp6 > A2B_B2A.bp6.rbh.pairs

-min_similarity       [0]
-max_lenDiffR         [0]
-log_lines            [0]

-gene1_list           [\$list_filename] A list with 1st column contains all required query gene IDs. 
-gene2_list           [\$list_filename] A list with 1st column contains all required subject gene IDs. 
-combine_gene12       [Boolean] Use both gene1_list and gene2_list as query/subject required IDs 
-paired_need          [Boolean] If this is given, -gene1_list and -gene2_list must have the same length. This won't work if geneX_list is empty. 

HH

defined $opts{'in_bp6'} or &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my %rbh_hash; 
$rbh_hash{'in_bp6'} = $opts{'in_bp6'}; 
$rbh_hash{'min_similarity'} = $opts{'min_similarity'}; 
$rbh_hash{'max_lenDiffR'} = $opts{'max_lenDiffR'}; 
$rbh_hash{'min_cscore'} = $opts{'min_cscore'}; 
$rbh_hash{'log_lines'}  = $opts{'log_lines'}; 
$rbh_hash{'gene1_need'} = &load_gene_list( $opts{'gene1_list'} ); 
$rbh_hash{'gene2_need'} = &load_gene_list( $opts{'gene2_list'} ); 
$rbh_hash{'paired_need'} = $opts{'paired_need'}; 
$rbh_hash{'combine_gene12'} = $opts{'combine_gene12'}; 

&SeqAlnSunhh::rbh_byBp6( \%rbh_hash ); 
my %used; 
&tsmsg("[Msg] Output RBH\n"); 
my %cnt; 
$cnt{'cntN_step'} = $opts{'log_lines'}; 
$cnt{'cnt'} = 0; 
for my $gA (sort keys %{$rbh_hash{'rbh'}}) {
	$cnt{'cnt'} ++; 
	$opts{'log_lines'} > 0 and &fileSunhh::log_section( $cnt{'cnt'}, \%cnt ) and &tsmsg("[Msg] Output $cnt{'cnt'} gene [$gA]\n"); 
	for my $gB (sort keys %{$rbh_hash{'rbh'}{$gA}}) {
		defined $used{$gA}{$gB} and next; 
		print STDOUT join("\t", $gA, $gB)."\n"; 
		$used{$gA}{$gB} = 1; 
		$used{$gB}{$gA} = 1; 
	}
}

&tsmsg("[Rec] $0 done.\n"); 


################################################################################
#    Subroutines 
################################################################################
sub load_gene_list {
	my ($fn) = @_; 
	my @back; 
	(defined $fn and $fn ne '') or return(\@back); 
	my $fh = &openFH($fn, '<'); 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		push(@back, $ta[0]); 
	}
	close($fh); 
	return(\@back); 
}# load_gene_list () 


