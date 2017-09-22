#!/usr/bin/perl
# Input must be sorted by position in order to save memory usage. 
# 2015-08-20 I want to write a script to filter sam reads for different purpose: 
#   Purpose 1 : I want uniquely-aligned read pairs with some (NM/read_length)% control; 
#   Purpose 2 : I want well-aligned read pairs with some (NM/read_length)% control; 
#                 Please note that here I don't filter insert size for well. 
#   Purpose 3 : I want best-aligned read pairs with some (NM/read_length)% control; 
#                 This is similar to -well_pair but filtering if XA_minNM < NM ; 
#   Purpose 4 : Both ends aligned, and filter the NM/XA_minNM ratio. 
#                 I don't care if the two ends map to the same chrID; 
#   Purpose 5 : Trim aligned reads in sam file. Trim the reads from their both ends to center by given length. 
#                 Please note thate in the new sam file, there will be no 'TAG:VALUE' pairs, and the mate-information is not correct ever. 
#   Purpose 6 : Filter sam alignments with global filtering rules. 
#   Purpose 7 : Add 'XT:i:##' information from in_picard_illuminaAdapterMarked.bam to input sam file, changing 'XT' to 'YT'; 
use strict; 
use warnings; 
use SeqAlnSunhh; 
use mathSunhh; 
# use IO::Handle; my $io = IO::Handle->new(); $io->fdopen(fileno(STDOUT), "w"); $io->autoflush(1);
use FileHandle; 
use IO::Handle; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	# Purpose class 
	'uniq_pair!', # Purpose 1. 
	'well_pair!', # Purpose 2. 
	'best_pair!', # Purpose 3. 
	'both_aln!',  # Purpose 4. 
	'trim_readEnds!', # Purpose 5. 
	'filter_sam!', # Purpose 6. 
	'add_XTi!', # Purpose 7. 
	
	# Local for 'well_pair'
	# Local for 'trim_readEnd'
	'trimLen:i', # Default 10 bp. 

	# Local for 'add_XTi'
	'bam_wiXTi:s', # The .bam file with 'XT:i:##' tag. 
	'tag4XTi:s',   # Default 'YT:i:', used to replace 'XT:i:' to avoid conflict with 'XT:A:?' ; 
	
	# Global filtering 
	'max_NM_ratio:f', # Default none. Recommend 0.01 for self_mapping. 
	'max_NM_ratio_inRdMatch:f', # Default none. Recommend 0.01 for self_mapping.
	'min_mapQ:i', # Default 0. 
	'noClip!', 
	'exe_samtools:s', # 
	'jar_picard:s', # 
	
	'help!', 
); 

STDOUT->autoflush(1); 

my %t_opts = %opts; 

# Global parameters. 
#my %flag_class; 
#$flag_class{'hDiff_Forward'} = &SeqAlnSunhh::mk_flag( 'keep'=>'0=1,2=0,3=0,4=0,5=1' , 'drop'=>'' ); 
#$flag_class{'hDiff_Reverse'} = &SeqAlnSunhh::mk_flag( 'keep'=>'0=1,2=0,3=0,4=1,5=0' , 'drop'=>'' ); 
$opts{'max_NM_ratio'} //= -1; 
$opts{'max_NM_ratio_inRdMatch'} //= -1; 
$opts{'min_mapQ'} //= 0; 

# Local parameters. 
$opts{'trimLen'} //= 10; 
$opts{'tag4XTi'} //= 'YT'; 


################################################################################
#      Invoke sub-routines. 
################################################################################
-t and !@ARGV and !( scalar(keys %t_opts) ) and &usage(); 

defined $opts{'help'} and &usage(); 

&tsmsg("[Rec] Begin $0\n"); 
if ( defined $opts{'uniq_pair'} or defined $opts{'well_pair'} or defined $opts{'best_pair'} ) {
	&get_uniq_pair(); 
} elsif ( defined $opts{'both_aln'} ) {
	&get_both_aln(); 
} elsif ( defined $opts{'trim_readEnds'} ) {
	&trim_readEnds(); 
} elsif ( defined $opts{'filter_sam'} ) {
	&filter_sam(); 
} elsif ( defined $opts{'add_XTi'} ) {
	&add_XTi(); 
}

&tsmsg("[Rec] Finish $0\n"); 

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
# -best_pair      [Boolean] Same to -well_pair but filtering if XA_minNM < NM ; 
#
# -both_aln       [Boolean] Both ends aligned, and filter the NM/XA_minNM ratio by parameter -max_NM_ratio ; 
#                   I don't care if the two ends map to the same chrID; 
#
# -trim_readEnds  [Boolean] Trim the reads alignment by removing the terminate bases from both ends. 
#   -trimLen      [$opts{'trimLen'}] The trimming length from each read end. 
#
# -filter_sam     [Boolean] Filter sam alignments with global filtering rules. 
#
# -add_XTi        [Boolean] Add 'XT:i:##' information for adapter information to input sam. 
#   -bam_wiXTi    [filename] Required with -add_XTi . The raw 'XT:i:##' exists here. 
#   -tag4XTi      [$opts{'tag4XTi'}] A string to replace 'XT:i:'. 
# 
# Global filtering: 
# -max_NM_ratio   [-1] Maximum NM/rdLen ratio accepted. 
# -max_NM_ratio_inRdMatch [-1] NM/Match_rdLen ratio. Only for filter_sam
# -min_mapQ       [0] Minimum mapping quality accepted. 
# -noClip         [Boolean] Do not allow hard/soft clip if given. 
# -exe_samtools   [String] for samtools path. 
# -jar_picard     [String] for picard.jar path. Preferred. 
# 
################################################################################
HH
	exit 1; 
}

sub add_XTi {
	defined $opts{'bam_wiXTi'} or &stopErr("[Err] -bam_wiXTi needed when using -add_XTi .\n"); 
	P1_in_add_XTi: 
	if (defined $opts{'jar_picard'}) {
		open(XTBAM, '-|', "java -jar $opts{'jar_picard'} ViewSam INPUT=$opts{'bam_wiXTi'} RECORDS_ONLY=true ") or &stopErr("[Err] $!\n"); 
	} elsif (defined $opts{'exe_samtools'}) {
		my $opt_S = ( $opts{'bam_wiXTi'} =~ m/\.sam$/i ) ? '-S' : ''; 
		open(XTBAM, '-|', "$opts{'exe_samtools'} view $opt_S $opts{'bam_wiXTi'}") or &stopErr("[Err] $!\n"); 
	} elsif (-e "$ENV{'HOME'}/bin/picard.jar") {
		$opts{'jar_picard'} = "$ENV{'HOME'}/bin/picard.jar"; 
		goto P1_in_add_XTi; 
	} else {
		&stopErr("[Err] At least need one of -jar_picard and -exe_samtools \n"); 
	}
	my %flag_R1 = map { $_ => 'R1' } keys %{ &SeqAlnSunhh::mk_flag( 'keep'=>'6=1' ) }; 
	my %flag_R2 = map { $_ => 'R2' } keys %{ &SeqAlnSunhh::mk_flag( 'keep'=>'7=1' ) }; 
	my %flag_Ru = map { $_ => 'Ru' } keys %{ &SeqAlnSunhh::mk_flag( 'keep'=>'6=0,7=0' ) }; 
	my %flag_R12 = (%flag_R1, %flag_R2, %flag_Ru); 
	my %rd2XTi; 
	while (<XTBAM>) {
		$. % 1e6 == 1 and &tsmsg("[Msg] Pre-reading $. reads.\n"); 
		chomp; 
		my @ta = split(/\t/, $_); 
		my $xt_tag = ''; 
		for (my $i=11; $i<@ta; $i++) {
			$ta[$i] =~ s!^XT:i:(\S+)$!$opts{'tag4XTi'}:i:$1!o or next; 
			$xt_tag = $ta[$i]; 
			last; 
		}
		$xt_tag eq '' and next; 
		defined $flag_R12{ $ta[1] } or &stopErr("[Err] Bad flag [$ta[1]]\n"); 
		$rd2XTi{$ta[0]}{ $flag_R12{$ta[1]} } = $xt_tag; 
	}
	close(XTBAM); 
	open(OO, '>', "/dev/stdout") or die "$!\n"; 
	select(OO); 
	$| = 1; 
	while (<>) {
		$. % 1e6 == 1 and &tsmsg("[Msg] Current sam $. line.\n"); 
		chomp; 
		my @ta = split(/\t/, $_); 
		scalar(@ta) >= 11 or do { $ta[0] =~ m/^@/ or &stopErr("[Err] Bad line : $_\n"); print OO "$_\n"; next; }; 
		if (defined $rd2XTi{$ta[0]}) {
			defined $flag_R12{$ta[1]} or &stopErr("[Err] Bad flag [$ta[1]]\n"); 
			my $r12 = $flag_R12{$ta[1]}; 
			defined $rd2XTi{$ta[0]}{$r12} and $_ .= "\t$rd2XTi{$ta[0]}{$r12}"; 
		}
		print OO "$_\n"; 
	}
	return (); 
}# sub add_XTi() 

sub filter_sam {
	while (<>) {
		$. % 1000e3 == 1 and &tsmsg("[Msg] Processing $. line.\n"); 
		chomp; 
		my @ta = split(/\t/, $_); 
		($ta[0] =~ m/^@/ and @ta <= 10) and do { print STDOUT "$_\n"; next; }; 
		my $sam_href = &SeqAlnSunhh::sam_line2hash(\@ta, []); 

		$sam_href->{'mapq'} >= $opts{'min_mapQ'} or do { next; }; 

		if ($opts{'max_NM_ratio_inRdMatch'} >= 0) {
			&SeqAlnSunhh::sam_hash_addKey($sam_href, ['NM', 'cigar_href']); 
			$sam_href->{'NM'} > $sam_href->{'cigar_href'}->{'MatchRdLen'} * $opts{'max_NM_ratio_inRdMatch'} and do { next; }; 
		}
		if ($opts{'max_NM_ratio'} >= 0) {
			&SeqAlnSunhh::sam_hash_addKey($sam_href, ['NM', 'read_len', 'cigar_href']); 
			$sam_href->{'NM'} + $sam_href->{'cigar_href'}{'Hlen'} + $sam_href->{'cigar_href'}{'Slen'} > $sam_href->{'read_len'} * $opts{'max_NM_ratio'} and do { next; }; 
		}
		if ($opts{'noClip'}) {
			&SeqAlnSunhh::sam_hash_addKey($sam_href, ['cigar_href']); 
			( $sam_href->{'cigar_href'}{'Hlen'} > 0 or $sam_href->{'cigar_href'}{'Slen'} > 0 ) and do { next; }; 
		}

		print STDOUT "$_\n"; 
	}
	return (); 
}# filter_sam() 

# Local  parameters: -uniq_pair
# Global parameters: -max_NM_ratio [:f] -min_mapQ [:i]
sub get_uniq_pair {
	my %bad_rd; 
	my @lines; 
	my $prev_chrID; 
	while (<>) {
		$. % 100e3 == 1 and &tsmsg("[Msg] Processing $. line.\n"); 
		# m/^\@/ and do { print; next; }; 
		
		chomp; 
		my @ta = split(/\t/, $_); 
		($ta[0] =~ m/^@/ and @ta <= 10) and do { print STDOUT "$_\n"; next; }; 
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
			$sam_href->{'NM'} > $sam_href->{'XA_minNM'} and do { $bad_rd{$rdID} = 1; next; }; 
		} elsif ( defined $opts{'well_pair'} ) {
		} else {
			&SeqAlnSunhh::sam_hash_addKey($sam_href, ['is_uniqBest']); 
			$sam_href->{'is_uniqBest'} == 1 or do { $bad_rd{$rdID} = 1; next; } ; 
		}
		if ($opts{'max_NM_ratio'} >= 0) {
			&SeqAlnSunhh::sam_hash_addKey($sam_href, ['NM', 'read_len', 'cigar_href']); 
			$sam_href->{'NM'} + $sam_href->{'cigar_href'}{'Hlen'} + $sam_href->{'cigar_href'}{'Slen'} > $sam_href->{'read_len'} * $opts{'max_NM_ratio'} and do { $bad_rd{$rdID} = 1; next; }; 
		}
		if ($opts{'noClip'}) {
			&SeqAlnSunhh::sam_hash_addKey($sam_href, ['cigar_href']); 
			( $sam_href->{'cigar_href'}{'Hlen'} > 0 or $sam_href->{'cigar_href'}{'Slen'} > 0 ) and do { $bad_rd{$rdID} = 1; next; }; 
		}
		
		push(@lines, [$rdID, $_]); 
	}
	if (@lines > 0) {
		&SeqAlnSunhh::print_sam_lines(\@lines, \%bad_rd); 
		@lines = (); 
		%bad_rd = (); 
	}
}# get_uniq_pair

sub get_both_aln {
	my %bad_rd; 
	my %rd_cnt_r1; 
	my %rd_cnt_r2; 
	my @lines; 
	my %paired_flag = %{ &SeqAlnSunhh::mk_flag('keep'=>'0=1,2=0,3=0') }; 
	my %r1_flag = %{ &SeqAlnSunhh::mk_flag('keep'=>'0=1,2=0,3=0,6=1,7=0') }; 
	my %r2_flag = %{ &SeqAlnSunhh::mk_flag('keep'=>'0=1,2=0,3=0,7=1,6=0') }; 
	while (<>) {
		$. % 100e3 == 1 and &tsmsg("[Msg] Processing $. line.\n"); 
		# m/^\@/ and do { print; next; }; 

		chomp; 
		my @ta = split(/\t/, $_); 
		($ta[0] =~ m/^@/ and @ta <= 10) and do { print STDOUT "$_\n"; next; }; 
		my $sam_href = &SeqAlnSunhh::sam_line2hash(\@ta ); 
		defined $paired_flag{ $ta[1] } or next; 
		my $rdID = $sam_href->{'qname'}; 
		defined $bad_rd{ $rdID } and next; 
		$sam_href->{'mapq'} >= $opts{'min_mapQ'} or do { $bad_rd{$rdID} = 1; next; }; 
		if ($opts{'max_NM_ratio'} >= 0) {
			&SeqAlnSunhh::sam_hash_addKey($sam_href, ['NM', 'read_len', 'XA_minNM', 'cigar_href']); 
			my $vv = &mathSunhh::min( $sam_href->{'NM'}, $sam_href->{'XA_minNM'} ); 
			$vv + $sam_href->{'cigar_href'}{'Hlen'} + $sam_href->{'cigar_href'}{'Slen'} > $sam_href->{'read_len'} * $opts{'max_NM_ratio'} and do { $bad_rd{$rdID} = 1; next; }; 
		}
		if ($opts{'noClip'}) {
			&SeqAlnSunhh::sam_hash_addKey($sam_href, ['cigar_href']); 
			( $sam_href->{'cigar_href'}{'Hlen'} > 0 or $sam_href->{'cigar_href'}{'Slen'} > 0 ) and do { $bad_rd{$rdID} = 1; next; }; 
		}
		

		if ( defined $r1_flag{$ta[1]} ) {
			$rd_cnt_r1{$rdID} ++; 
		} elsif ( defined $r2_flag{$ta[1]} ) {
			$rd_cnt_r2{$rdID} ++; 
		} else {
			&stopErr("[Err] Bad flag [$ta[1]]\n"); 
		}
		push(@lines, [$rdID, $_]); 
	}
	&tsmsg("[Msg] Filtering reads for pairing\n"); 
	if (@lines > 0) {
		for my $rdID (keys %rd_cnt_r1) {
			defined $bad_rd{$rdID} and next; 
			defined $rd_cnt_r2{$rdID} or $bad_rd{$rdID} = 1; 
		}
		&SeqAlnSunhh::print_sam_lines(\@lines, \%bad_rd); 
		@lines = (); 
		%bad_rd = (); 
		%rd_cnt_r1 = %rd_cnt_r2 = (); 
	}
}

sub trim_readEnds {
	my %flag_fwd = %{ &SeqAlnSunhh::mk_flag('keep'=>'2=0,4=0') }; 
	my %flag_rev = %{ &SeqAlnSunhh::mk_flag('keep'=>'2=0,4=1') }; 
	while (<>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		($ta[0] =~ m/^@/ and @ta <= 10) and do { print STDOUT "$_\n"; next; }; 
		$ta[4] >= $opts{'min_mapQ'} or next; 
		if ( $opts{'max_NM_ratio'} >= 0 ) {
			my $sam_href = &SeqAlnSunhh::sam_line2hash( \@ta, ['NM', 'read_len', 'XA_minNM', 'cigar_href'] ); 
			my $vv = &mathSunhh::min( $sam_href->{'NM'}, $sam_href->{'XA_minNM'} ); 
			$vv + $sam_href->{'cigar_href'}{'Hlen'} + $sam_href->{'cigar_href'}{'Slen'} > $sam_href->{'read_len'} * $opts{'max_NM_ratio'} and next; 
		}
		(defined $flag_fwd{$ta[1]} or defined $flag_rev{$ta[1]}) or do { print STDOUT join("\t", @ta[0 .. 10])."\n"; next; }; 
		my ($rd_len, $ref_span) = &SeqAlnSunhh::cigar_array2len( &SeqAlnSunhh::cigar_str2array( $ta[5] ) );
		my ($new_cigar_aref, $new_drop_lp, $new_drop_rp, $left_mv_refP, $right_mv_refP) = &SeqAlnSunhh::trim_cigar_str_bothEnd($ta[5], $opts{'trimLen'}); 
		scalar( @$new_cigar_aref ) > 0 or next; 
		my $new_cigar_str = &SeqAlnSunhh::cigar_array2str( $new_cigar_aref ); 
		my $new_rd_bp = substr( $ta[9],  $new_drop_lp, $rd_len-$new_drop_lp-$new_drop_rp ); 
		my $new_rd_qb = substr( $ta[10], $new_drop_lp, $rd_len-$new_drop_lp-$new_drop_rp ); 
		my $new_refP = $ta[3]+$left_mv_refP; 
		print STDOUT join("\t", @ta[0 .. 2], $new_refP, $ta[4], $new_cigar_str, @ta[6 .. 8], $new_rd_bp, $new_rd_qb)."\n"; 
	}
}# sub trim_readEnds 

################################################################################
#      Inner Sub-routines
################################################################################

