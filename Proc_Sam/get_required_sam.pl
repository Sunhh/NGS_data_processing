#!/usr/bin/perl
# Input must be sorted by position in order to save memory usage. 
# 2015-08-20 I want to write a script to filter sam reads for different purpose: 
#   Purpose 1 : I want uniquely-aligned read pairs with some (NM/read_length)% control; 
#   Purpose 2 : I want well-aligned read pairs with some (NM/read_length)% control; 
#                 Please note that here I don't filter insert size for well. 
#   Purpose 3 : I want best-aligned read pairs with some (NM/read_length)% control; 
#                 This is similar to -well_pair but filtering if XA_minNM <= NM ; 
use strict; 
use warnings; 
use SeqAlnSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	# Purpose class 
	'uniq_pair!', # Purpose 1. 
	'well_pair!', # Purpose 2. 
	'best_pair!', # Purpose 3. 
	
	# Local for 'well_pair'
	
	# Global filtering 
	'max_NM_ratio:f', # Default none. Recommend 0.01 for self_mapping. 
	'min_mapQ:i', # Default 0. 
	
	'help!', 
); 

-t and !@ARGV and !( scalar(keys %opts) ) and &usage(); 

# Global parameters. 
#my %flag_class; 
#$flag_class{'hDiff_Forward'} = &SeqAlnSunhh::mk_flag( 'keep'=>'0=1,2=0,3=0,4=0,5=1' , 'drop'=>'' ); 
#$flag_class{'hDiff_Reverse'} = &SeqAlnSunhh::mk_flag( 'keep'=>'0=1,2=0,3=0,4=1,5=0' , 'drop'=>'' ); 
$opts{'max_NM_ratio'} //= -1; 
$opts{'min_mapQ'} //= 0; 


################################################################################
#      Invoke sub-routines. 
################################################################################
defined $opts{'help'} and &usage(); 
if ( defined $opts{'uniq_pair'} or defined $opts{'well_pair'} or defined $opts{'best_pair'} ) {
	&get_uniq_pair(); 
}


################################################################################
#      Sub-routines to execute 
################################################################################
sub usage {
	print STDOUT <<HH; 
################################################################################
# perl $0 in_pair_aln.sam 
# 
# -help 
# 
# Stratergy: 
# -uniq_pair      [Boolean] read pairs with both ends uniquely-aligned to the same chr ; 
#                   Applies: -max_NM_ratio -min_mapQ 
# 
# -well_pair      [Boolean] well-aligned read pairs with some (NM/read_length)% control; 
# -best_pair      [Boolean] Same to -well_pair but filtering if XA_minNM <= NM ; 
# 
# Global filtering: 
# -max_NM_ratio   [-1] Maximum NM/rdLen ratio accepted. 
# -min_mapQ       [0] Minimum mapping quality accepted. 
# 
################################################################################
HH
	exit 1; 
}

# Local  parameters: -uniq_pair
# Global parameters: -max_NM_ratio [:f] -min_mapQ [:i]
sub get_uniq_pair {
	my %bad_rd; 
	my @lines; 
	my $prev_chrID; 
	while (<>) {
		$. % 100e3 == 1 and &tsmsg("[Msg] Processing $. line.\n"); 
		m/^\@/ and do { print; next; }; 
		
		chomp; 
		my @ta = split(/\t/, $_); 
		# This step takes 3 times duration compared to the original one because of the TAG-to-hash calculation. 
		# But this time increase is acceptable because the script is easier to manage. 
		my $sam_href = &SeqAlnSunhh::sam_line2hash(\@ta, ['hDiff_Pair']); 
		my $rdID = $sam_href->{'qname'}; 
		defined $bad_rd{ $rdID } and next; 
		$sam_href->{'hDiff_Pair'} or next; 
		$sam_href->{'rnext'} eq '=' or next; 
		$sam_href->{'mapq'} >= $opts{'min_mapQ'} or do { $bad_rd{$rdID} = 1; next; }; 
		$prev_chrID //= $sam_href->{'rname'}; 
		if ( $prev_chrID ne $sam_href->{'rname'} ) {
			&SeqAlnSunhh::print_sam_lines(\@lines, \%bad_rd); 
			@lines=(); 
			%bad_rd=(); 
			$prev_chrID = $sam_href->{'rname'}; 
		}
		if ( defined $opts{'best_pair'} ) {
			&SeqAlnSunhh::sam_hash_addKey($sam_href, ['NM', 'XA_minNM']); 
			$sam_href->{'NM'} >= $sam_href->{'XA_minNM'} and do { $bad_rd{$rdID} = 1; next; }; 
		} elsif ( defined $opts{'well_pair'} ) {
		} else {
			&SeqAlnSunhh::sam_hash_addKey($sam_href, ['is_uniqBest']); 
			$sam_href->{'is_uniqBest'} == 1 or do { $bad_rd{$rdID} = 1; next; } ; 
		}
		if ($opts{'max_NM_ratio'} >= 0) {
			&SeqAlnSunhh::sam_hash_addKey($sam_href, ['NM', 'read_len']); 
			$sam_href->{'NM'} > $sam_href->{'read_len'} * $opts{'max_NM_ratio'} and do { $bad_rd{$rdID} = 1; next; }; 
		}
		
		push(@lines, [$rdID, $_]); 
	}
	if (@lines > 0) {
		&SeqAlnSunhh::print_sam_lines(\@lines, \%bad_rd); 
		@lines = (); 
		%bad_rd = (); 
	}
}# get_uniq_pair

################################################################################
#      Inner Sub-routines
################################################################################

