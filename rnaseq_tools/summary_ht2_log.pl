#!/usr/bin/perl
use strict; 
use warnings; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"getLibInfo!", 
); 

!@ARGV and die "perl $0 o_prefix proc_rnaseq/LOG/Pub3_LFMOCK_Rep1_time2.ht2.log > proc_rnaseq/LOG/Pub3_LFMOCK_Rep1_time2.ht2.log.tbl\n"; 

my $opref = shift; 

my %h; 
while (<>) {
	chomp; 
	if      (m!^(\d+) reads; of these!) {
		$h{'totRd'} = $1; 
	} elsif (m!^\s+(\d+) \(\S+ were paired!) {
		$h{'type'} = 'PE'; 
	} elsif (m!^\s+(\d+) \(\S+ were unpaired!) {
		$h{'type'} = 'SE'; 
	} elsif (m!^\s+(\d+) \(\S+ aligned concordantly 0 times!) {
		; 
	} elsif (m!^\s+(\d+) \(\S+ aligned concordantly exactly 1 time!) {
		$h{'pe_goodUmap'} = $1; 
	} elsif (m!^\s+(\d+) \(\S+ aligned concordantly \>1 times!) {
		$h{'pe_goodMmap'} = $1; 
	} elsif (m!^\s+\-+\s*$!) {
		; 
	} elsif (m!^\s+(\d+) pairs aligned concordantly 0 times!) {
		; 
	} elsif (m!^\s+(\d+) \(\S+ aligned discordantly 1 time!) {
		$h{'pe_discordMap'} = $1; 
	} elsif (m!^\s+(\d+) pairs aligned 0 times concordantly or discordantly!) {
		; 
	} elsif (m!^\s+(\d+) mates make up the pairs!) {
		$h{'pe_notP'} = $1; 
	} elsif (m!^\s+(\d+) \(\S+ aligned 0 times!) {
		; 
	} elsif (m!^\s+(\d+) \(\S+ aligned exactly 1 time!) {
		$h{'pe_notP_isS_U'} = $1; 
		$h{'se_isS_U'} = $1; 
	} elsif (m!^\s+(\d+) \(\S+ aligned \>1 times!) {
		$h{'pe_notP_isS_M'} = $1; 
		$h{'se_isS_M'} = $1; 
	} elsif (m!^([\d.]+\% overall alignment rate)!) {
		; 
	} else {
		die "unknown line: $_\n"; 
	}
}

for (qw/totRd pe_goodUmap pe_goodMmap pe_discordMap pe_notP pe_notP_isS_U se_isS_U se_isS_U pe_notP_isS_M se_isS_M/) {
	$h{$_} //= 0; 
}
defined $h{'type'} or die "unknown type\n"; 
if ( $h{'type'} eq 'PE' ) {
	$h{'totRd'} *= 2; 
	$h{'pe_goodUmap'} *= 2; 
	$h{'pe_goodMmap'} *= 2; 
	$h{'pe_discordMap'} *= 2; 
	$h{'pe_notP_isS_U'} *= 2; 
	$h{'pe_notP_isS_M'} *= 2; 
}
$h{'AllMapRd'} = &get_sum( @h{qw/pe_goodUmap pe_goodMmap pe_discordMap pe_notP_isS_U pe_notP_isS_M/} ); 
$h{'UMapRd'}   = &get_sum(@h{qw/pe_goodUmap pe_discordMap pe_notP_isS_U/}); 

if ($opts{'getLibInfo'}) {
	# Could have value_type : 
	#   lib_rd_len        : Read length ; 
	#   lib_tNSum_woRdLen : -1 indicates I need to count ; 
	#   lib_sumMapRd      : Should be background total fragments number; 
	# $h{'type'} eq 'PE' and do { $h{'AllMapRd'} /= 2; $h{'totRd'} /= 2;  }; 
	print STDOUT join("\t", '#libID', qw/value_type      value/)."\n"; 
	print STDOUT join("\t", $opref,   'lib_seqType',     $h{'type'})."\n"; 
	if ($h{'type'} eq 'PE') {
		print STDOUT join("\t", $opref,   'lib_sumMapRd',    $h{'AllMapRd'}/2 )."\n"; 
		print STDOUT join("\t", $opref,   'lib_totalRdP',    $h{'totRd'}/2 )."\n"; 
	} elsif ($h{'type'} eq 'SE') {
		print STDOUT join("\t", $opref,   'lib_sumMapRd',    $h{'AllMapRd'} )."\n"; 
		print STDOUT join("\t", $opref,   'lib_totalRdP',    $h{'totRd'} )."\n"; 
	} else {
		die "Bad lib_seqType [$h{'type'}]\n"; 
	}
} else {
	print STDOUT join("\t", qw/LibID SeqType Total_RdN PE_concord_U PE_concord_M PE_discord SE_U SE_M AllMapRd Perc_AllMapRat UMapRd Perc_UMapRat/)."\n"; 
	print STDOUT join("\t", 
		$opref, 
		$h{'type'},  # SeqType
		$h{'totRd'}, # Total_RdN 
		$h{'pe_goodUmap'}, # PE_concord_U
		$h{'pe_goodMmap'}, # PE_concord_M
		$h{'pe_discordMap'}, # PE_discord
		$h{'pe_notP_isS_U'}, # SE_U 
		$h{'pe_notP_isS_M'}, # SE_M 
		$h{'AllMapRd'}, 
		sprintf("%.2f", 100 * $h{'AllMapRd'} / $h{'totRd'} ),  # Perc_AllMapRat = (PE_concord_U + PE_concord_M + PE_discord + SE_U + SE_M)/Total_RdN
		$h{'UMapRd'}, 
		sprintf("%.2f", 100 * $h{'UMapRd'}   / $h{'totRd'} ) # Perc_UMapRat   = (PE_concord_U + PE_discord + SE_U )/Total_RdN
	)."\n"; 
}

sub get_sum {
	my $v=0; 
	for (@_) {
		$v+=$_; 
	}
	return($v); 
}

# [Sunhh@bioinfor01 03_alignToWM97ngsV1Scf_wiHiSat2]$ less -S proc_rnaseq/LOG/Pub3_LFMOCK_Rep1_time1.ht2.log
# 16091213 reads; of these:
#   16091213 (100.00%) were paired; of these:
#     677557 (4.21%) aligned concordantly 0 times
#     15010551 (93.28%) aligned concordantly exactly 1 time
#     403105 (2.51%) aligned concordantly >1 times
#     ----
#     677557 pairs aligned concordantly 0 times; of these:
#       142590 (21.04%) aligned discordantly 1 time
#     ----
#     534967 pairs aligned 0 times concordantly or discordantly; of these:
#       1069934 mates make up the pairs; of these:
#         636808 (59.52%) aligned 0 times
#         414517 (38.74%) aligned exactly 1 time
#         18609 (1.74%) aligned >1 times
# 98.02% overall alignment rate
# 
# 
# [Sunhh@bioinfor01 03_alignToWM97ngsV1Scf_wiHiSat2]$ less -S proc_rnaseq/LOG/self_wpW_LF10_a_Rep2_time1.ht2.log
# 8851215 reads; of these:
#   8851215 (100.00%) were unpaired; of these:
#     152788 (1.73%) aligned 0 times
#     8553288 (96.63%) aligned exactly 1 time
#     145139 (1.64%) aligned >1 times
# 98.27% overall alignment rate



