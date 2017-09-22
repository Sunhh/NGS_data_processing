#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"rawScfKeyFile:s", 
); 
# [Sunhh@whale super_scaffold]$ head -3 SP_Nov_11_13_30_8_90_3_superscaffold.fasta_contig.bed
# 1       1       126451  Contig_1
# 1       127456  180114  Contig_2
# 1       180182  209486  Contig_3
#
my $help_txt = <<HH; 

perl $0 In_superscaffold_merged_BspQI_key.txt In_superscaffold_merged.agp > In_superscaffold_merged.bed

-help 

-rawScfKeyFile      [''] In_raw_scf_BspQI_key.txt

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
!@ARGV and &LogInforSunhh::usage($help_txt); 

my %key_info = %{ &load_key_txt( $ARGV[0] ) }; 
my @txt_arr = @{ &load_agp_file( $ARGV[1] ) }; 
my %key_info_raw; 
defined $opts{'rawScfKeyFile'} and %key_info_raw = %{ &load_key_txt( $opts{'rawScfKeyFile'} ) }; 

for my $la (@txt_arr) {
	$la->[4] =~ m/^(N|U)$/ and next; 
	my $compID = $key_info{'Raw2Comp'}{$la->[0]}; 
	my $bed_eleID = $la->[5]; 
	defined $opts{'rawScfKeyFile'} and defined $key_info_raw{'Raw2Comp'}{$la->[5]} and $bed_eleID="$key_info_raw{'Raw2Comp'}{$la->[5]}:$bed_eleID"; 
	print join("\t", $compID, $la->[1], $la->[2], $bed_eleID)."\n"; 
}


# [Sunhh@whale super_scaffold]$ head -3 SP_Nov_11_13_30_8_90_3_superscaffold_merged.agp
# ##agp-version   2.0
# Super_scaffold_248      1       341379  1       W       SpoScf_00639    1       341379  -
# Super_scaffold_248      341380  600696  2       N       259317  scaffold        yes     map
sub load_agp_file {
	# https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/ 
	my $fh = &openFH($_[0], '<'); 
	my @back; 
	while (<$fh>) {
		m/^\s*(#|$)/ and next; 
		chomp; 
		my @ta = split(/\t/, $_); 
		push(@back, [@ta]); 
	}
	close($fh); 
	return(\@back); 
}

# [Sunhh@whale super_scaffold]$ head -5 SP_Nov_11_13_30_8_90_3_superscaffold_BspQI_key.txt
# CMAP = /data/Sunhh/spinach/16.bionano/01.stitch/Nov_11_r2/default_alignment/stitch3/SP_Nov_11_13_30_8_90_3_superscaffold_BspQI.cmap
# # filter: Minimum Labels = 5
# # filter: Minimum Size (Kb) = 20
# CompntId        CompntName      CompntLength
# 1       Super_scaffold_248      3853603
sub load_key_txt {
	my $fh = &openFH($_[0], '<'); 
	my %h; 
	while (<$fh>) {
		m/^\s*(#|$)/ and next; 
		chomp; 
		my @ta = split(/\t/, $_); 
		$ta[0] eq 'CompntId' and next; 
		$h{'Comp2Raw'}{$ta[0]} = $ta[1]; 
		$h{'Raw2Comp'}{$ta[1]} = $ta[0]; 
		if (defined $ta[2]) {
			$h{'Raw2Len'}{$ta[1]} = $ta[2]; 
			$h{'Comp2Len'}{$ta[0]} = $ta[2]; 
		}
	}
	close($fh); 
	return (\%h); 
}

