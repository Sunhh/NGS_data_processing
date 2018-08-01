#!/usr/bin/perl
# 2018-07-30 : Fixing paired-end read pairs' alignments; 
use strict; 
use warnings; 
use LogInforSunhh; 
use SeqAlnSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"inBam:s", # 
	"outBam:s", 
	"nSorted!", 
	"para_samTsrt:s", 
	"exe_samtools:s", 
	"fixSRRPairTag!", 
	"help!", 
); 

$opts{'exe_samtools'} //= 'samtools'; 
$opts{'para_samTsrt'} //= ' -@ 4 -m 5G '; 

my %str2flag; 

my $help_txt = <<HH; 
################################################################################
# perl $0   -inBam hisat2_filtered.bam -outBam hisat2_fixNH.bam 
#
# -help 
#
# -exe_samtools   ['$opts{'exe_samtools'}']
# -para_samTsrt   ['$opts{'para_samTsrt'}']
# -nSorted        [Boolean] Already sorted by name if given. 
#
# -fixSRRPairTag  [Boolean] If given, I suppose the input paired reads IDs are *.1 / *.2; 
#                           And I will remove .1/.2 information; 
################################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
defined $opts{'inBam'}  or &LogInforSunhh::usage($help_txt); 
defined $opts{'outBam'} or &LogInforSunhh::usage($help_txt); 


my %flag_RS = %{ &SeqAlnSunhh::mk_flag( 'keep'=>'0=0,6=0,7=0' ) }; 
my %flag_R1 = %{ &SeqAlnSunhh::mk_flag( 'keep'=>'0=1,6=1,7=0' ) }; 
my %flag_R2 = %{ &SeqAlnSunhh::mk_flag( 'keep'=>'0=1,6=0,7=1' ) }; 
my %flag_RP = %{ &SeqAlnSunhh::mk_flag( 'keep'=>'0=1' ) };
my %flag_mateMap = %{ &SeqAlnSunhh::mk_flag( 'keep'=>'3=0' ) }; 
my %flag_readMap = %{ &SeqAlnSunhh::mk_flag( 'keep'=>'2=0' ) }; 
my %flag2R; 
for (keys %flag_RS) { $flag2R{$_} = 'RS'; }
for (keys %flag_R1) { $flag2R{$_} = 'R1'; }
for (keys %flag_R2) { $flag2R{$_} = 'R2'; }
my %flag_supp    = %{ &SeqAlnSunhh::mk_flag( 'keep'=>'11=1' ) }; 
my %flag_notSupp = %{ &SeqAlnSunhh::mk_flag( 'keep'=>'11=0' ) }; 
my %flag2supp; 
for (keys %flag_supp    ) { $flag2supp{$_} = '1'; }
for (keys %flag_notSupp ) { $flag2supp{$_} = '0'; }
my %flag_primary = %{ &SeqAlnSunhh::mk_flag( 'keep'=>'8=0' ) }; 
my %flag_notPrim = %{ &SeqAlnSunhh::mk_flag( 'keep'=>'8=1' ) }; 
my %flag2primary; 
for (keys %flag_primary) { $flag2primary{$_} = '1'; }
for (keys %flag_notPrim) { $flag2primary{$_} = '0'; }

my $tmpDir = &fileSunhh::new_tmp_dir('create'=>1); 

my $baseFn = &fileSunhh::_basename( $opts{'inBam'} ); 
if ( $opts{'nSorted'} ) {
	open F,'-|', "$opts{'exe_samtools'} view -h $opts{'inBam'}" or die "$!\n"; 
	&fileSunhh::_rmtree($tmpDir); 
} else {
	&exeCmd_1cmd("$opts{'exe_samtools'} sort -n $opts{'para_samTsrt'} -o $tmpDir/nSrt_$baseFn $opts{'inBam'}"); 
	open F,'-|', "$opts{'exe_samtools'} view -h $tmpDir/nSrt_$baseFn" or die "$!\n"; 
}
open O,'|-', "$opts{'exe_samtools'} view -o $opts{'outBam'} -" or die "$!\n"; 
my @sameRdPS = ('', '', {}); # ("readID\tread_RP/S", "P/S", {  })
# 2018-07-30 : Sometimes, one end of a pair will be eleminate, so its mate-pair's mate-alignment needs to be changed to unmapped. 
#   So &sam_flag_infor needs to be used. 
my %tmp_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>1e6 ); 
while (<F>) {
	&fileSunhh::log_section($., \%tmp_cnt) and &tsmsg("[Msg] Processing $. line in [$opts{'inBam'}].\n"); 
	if (m!^\@!) {
		print O; 
		next; 
	}
	chomp; 
	my @ta=split(/\t/, $_); 
	my $rS12 = $flag2R{$ta[1]}; 
	( defined $rS12 and $rS12 ne '' ) or die "Bad flag [$ta[1]]\n"; 
	if ($opts{'fixSRRPairTag'}) {
		if ( $rS12 eq 'R1' ) {
			$ta[0] =~ s!\.1$!! or &stopErr("[Err] Failed to remove SRR tag .1 from read name [$ta[0]]\n"); 
		} elsif ( $rS12 eq 'R2' ) {
			$ta[0] =~ s!\.2$!! or &stopErr("[Err] Failed to remove SRR tag .2 from read name [$ta[0]]\n"); 
		} elsif ( $rS12 eq 'RS' ) {
			$ta[0] =~ s!\.1$!! or &stopErr("[Err] Failed to remove SRR tag .1 from read name [$ta[0]]\n"); 
		} else {
			&stopErr("[Err] Unkown RS12 [$rS12]\n"); 
		}
	}
	my $rSP  = ($rS12 eq 'RS') ? 'S' : 'P' ; 
	my $rdKey = "$ta[0]\t$rSP"; 
	for (my $i=11; $i<@ta; $i++ ) {
		$ta[$i] =~ m!^NH:i:\d+$! or next; 
		splice(@ta, $i, 1); 
		last; 
	}
	if ($sameRdPS[0] eq '') {
		@sameRdPS = ( $rdKey, $rSP, {} ); 
		push(@{$sameRdPS[2]{$rS12}{'line'}}, [@ta]); 
		defined $flag_readMap{$ta[1]} and $ta[2] ne '*' and $sameRdPS[2]{$rS12}{'NH'} ++; 
	} elsif ( $sameRdPS[0] eq $rdKey ) {
		push(@{$sameRdPS[2]{$rS12}{'line'}}, [@ta]); 
		defined $flag_readMap{$ta[1]} and $ta[2] ne '*' and $sameRdPS[2]{$rS12}{'NH'} ++; 
	} else {
		&outSameRd(\@sameRdPS); 
		@sameRdPS = ( $rdKey, $rSP, {} ); 
		push(@{$sameRdPS[2]{$rS12}{'line'}}, [@ta]); 
		defined $flag_readMap{$ta[1]} and $ta[2] ne '*' and $sameRdPS[2]{$rS12}{'NH'} ++; 
	}
}
&outSameRd(\@sameRdPS); 
@sameRdPS = (); 
close O; 
close F; 

$opts{'nSorted'} or &fileSunhh::_rmtree($tmpDir); 

&tsmsg("[Rec] Finished $opts{'inBam'}\n"); 

sub outSameRd {
	my ($ar) = @_; 
	my @sameRdPS = @$ar; 
	if ($sameRdPS[1] eq 'S') {
		for my $tline (@{$sameRdPS[2]{'RS'}{'line'}}) {
			if (defined $flag_readMap{$tline->[1]}) {
				print O join("\t", @$tline, "NH:i:$sameRdPS[2]{'RS'}{'NH'}")."\n"; 
			} else {
				print O join("\t", @$tline)."\n"; 
			}
		}
	} elsif ($sameRdPS[1] eq 'P') {
		my $nh = 0; 
		my %used_R; 

		if (defined $sameRdPS[2]{'R1'}) {
			for (my $i1=0; $i1<@{$sameRdPS[2]{'R1'}{'line'}}; $i1++) {
				my $tline_1 = $sameRdPS[2]{'R1'}{'line'}[$i1]; 
				if (defined $flag_readMap{$tline_1->[1]}) {
				} else {
					if (defined $flag_mateMap{$tline_1->[1]}) {
						# Check if there is a good R2 alignment; 
						my $R2_chrID_byR1 = ($tline_1->[6] eq '=') ? $tline_1->[2] : $tline_1->[6] ; 
						my $has_R2 = 0; 
						if (defined $sameRdPS[2]{'R2'}) {
							for (my $j2=0; $j2<@{$sameRdPS[2]{'R2'}{'line'}}; $j2++) {
								my $tline_2 = $sameRdPS[2]{'R2'}{'line'}[$j2]; 
								defined $flag_readMap{$tline_2->[1]} or  next; 
								defined $flag_mateMap{$tline_2->[1]} and next; 
								($tline_2->[2] eq $R2_chrID_byR1 and $tline_2->[3] eq $tline_1->[7]) or next; 
								defined $used_R{'R2'}{$j2} or $nh++; 
								$used_R{'R2'}{$j2} = 1; 
								$has_R2 = 1; 
							}
						}
						if ($has_R2 == 0) {
							# The mate is changed to unmapped. 
							my $R1_flag_ar = &SeqAlnSunhh::sam_flag_infor( $tline_1->[1] ); 
							$R1_flag_ar->[1][0] = 0; # the read is mapped in a proper pair
							$R1_flag_ar->[3][0] = 1; # the mate is unmapped
							$R1_flag_ar->[5][0] = 0; # strand of the mate (1 for reverse)
							$tline_1->[1] = &get_flag( $R1_flag_ar ); 
							$tline_1->[2] = '*'; 
							$tline_1->[3] = '0'; 
							$tline_1->[4] = '0'; 
							$tline_1->[5] = '*'; 
							$tline_1->[6] = '*'; 
							$tline_1->[7] = '0'; 
							$tline_1->[8] = '0'; 
						}
					}# if mate is mapped. 
					next; 
				}
				$tline_1->[2] eq '*' and next; 
				# If R1 is mapped, mapping to one position counts $nh only once. 
				defined $used_R{'R1'}{$i1} or $nh++; $used_R{'R1'}{$i1} = 1; 
				if (defined $flag_mateMap{$tline_1->[1]} and $tline_1->[6] ne '*' and $tline_1->[7] ne '*') {
					my $R2_chrID_byR1 = ($tline_1->[6] eq '=') ? $tline_1->[2] : $tline_1->[6] ; 
					my $has_R2 = 0; 
					if (defined $sameRdPS[2]{'R2'}) {
						for (my $j2=0; $j2<@{$sameRdPS[2]{'R2'}{'line'}}; $j2++) {
							defined $used_R{'R2'}{$j2} and next; # 1-to-1 relationship; 
							my $tline_2 = $sameRdPS[2]{'R2'}{'line'}[$j2]; 
							defined $flag_readMap{$tline_2->[1]} or next; 
							defined $flag_mateMap{$tline_2->[1]} or next; 
							# To look for a mate alignment, we requires : 
							#   |R1_ins|   == |R2_ins| [if both aligned]
							#   Both are supplementary or not; 
							#   Both are primary or not; 
							#   R2_ChrID == R2_ChrID_fromR1mate 
							#   R2_Pos   == R2_Pos_fromR1mate
							#   R1_Pos   == R1_Pos_fromR2mate 
							#   R1_ChrID == R1_ChriD_fromR2mate
							abs($tline_2->[8]) == abs($tline_1->[8]) or next; 
							$flag2supp{ $tline_2->[1] } eq $flag2supp{ $tline_1->[1] } or next; 
							$flag2primary{ $tline_2->[1] } eq $flag2primary{ $tline_1->[1] } or next; 
							( $tline_2->[2] eq $R2_chrID_byR1 and $tline_2->[3] eq $tline_1->[7] ) or next; 
							my $R1_chrID_byR2 = ( $tline_2->[6] eq '=') ? $tline_2->[2] : $tline_2->[6] ; 
							( $tline_1->[2] eq $R1_chrID_byR2 and $tline_1->[3] eq $tline_2->[7] ) or next; 
							# The mate R2 alignment should not count for $nh any more. 
							$used_R{'R2'}{$j2} = 1; 
							$has_R2 = 1; # There is a mate alignment. 
							# I don't exit at once because there are sometimes exactly same supplementary alignments in R2. 
							#   All of them need to be tagged as good_used; 
							last; # Using index; 
						}
					}
					if ($has_R2 == 0) {
						# The mate is changed to unmapped; 
						my $R1_flag_ar = &SeqAlnSunhh::sam_flag_infor( $tline_1->[1] ); 
						$R1_flag_ar->[1][0] = 0; # the read is mapped in a proper pair
						$R1_flag_ar->[3][0] = 1; # the mate is unmapped
						$R1_flag_ar->[5][0] = 0; # strand of the mate (1 for reverse)
						$tline_1->[1] = &get_flag( $R1_flag_ar ); 
						$tline_1->[6] = '='; 
						$tline_1->[7] = $tline_1->[3]; 
						$tline_1->[8] = '0'; 
					}
				}
			}
		}
		if (defined $sameRdPS[2]{'R2'}) {
			for (my $j2=0; $j2<@{$sameRdPS[2]{'R2'}{'line'}}; $j2++) {
				my $tline_2 = $sameRdPS[2]{'R2'}{'line'}[$j2]; 
				if ( defined $flag_readMap{$tline_2->[1]} ) {
				} else {
					if (defined $flag_mateMap{$tline_2->[1]}) {
						# Check if R1 alignment is a good mate. 
						defined $used_R{'R2'}{$j2} and next; 
						# The mate infor should be changed to unmapped; 
						$used_R{'R2'}{$j2} = 1; 

						my $R2_flag_ar = &SeqAlnSunhh::sam_flag_infor( $tline_2->[1] ); 
						$R2_flag_ar->[1][0] = 0; 
						$R2_flag_ar->[3][0] = 1; 
						$R2_flag_ar->[5][0] = 0; 
						$tline_2->[1] = &get_flag( $R2_flag_ar ); 
						$tline_2->[2] = '*'; 
						$tline_2->[3] = '0'; 
						$tline_2->[4] = '0'; 
						$tline_2->[5] = '*'; 
						$tline_2->[6] = '*'; 
						$tline_2->[7] = '0'; 
						$tline_2->[8] = '0'; 
					}
					next; 
				}
				$tline_2->[2] eq '*' and next; 
				my $tk_R2 = join("\t", @{$tline_2}[1,2,3,5,6,7,8]); 
				if (defined $used_R{'R2'}{$j2}) {
					# If used_R, this alignment has been counted. 
				} else {
					# This R2 alignment has not been counted, which means it lacks R1 alignment. 
					$nh++; 
					$used_R{'R2'}{$j2} = 1; 
					if (defined $flag_mateMap{$tline_2->[1]} and $tline_2->[6] ne '*' and $tline_2->[7] ne '*') {
						# The mate (R1) alignment has been removed. 
						my $R2_flag_ar = &SeqAlnSunhh::sam_flag_infor( $tline_2->[1] ); 
						$R2_flag_ar->[1][0] = 0; 
						$R2_flag_ar->[3][0] = 1; 
						$R2_flag_ar->[5][0] = 0; 
						$tline_2->[1] = &get_flag( $R2_flag_ar ); 
						$tline_2->[6] = '='; 
						$tline_2->[7] = $tline_2->[3]; 
						$tline_2->[8] = '0'; 
					}
				}
			}
		}
		
		if ( defined $sameRdPS[2]{'R1'}{'line'} ) {
			for my $tline_1 (@{$sameRdPS[2]{'R1'}{'line'}}) {
				if ( defined $flag_readMap{$tline_1->[1]} ) {
					print O join("\t", @$tline_1, "NH:i:$nh")."\n"; 
				} else {
					print O join("\t", @$tline_1)."\n"; 
				}
			}
		}
		if ( defined $sameRdPS[2]{'R2'}{'line'} ) {
			for my $tline_2 (@{$sameRdPS[2]{'R2'}{'line'}}) {
				if ( defined $flag_readMap{$tline_2->[1]} ) {
					print O join("\t", @$tline_2, "NH:i:$nh")."\n"; 
				} else {
					print O join("\t", @$tline_2)."\n"; 
				}
			}
		}
	} else {
		&stopErr("[Err] Bad read RP tag [$ar->[1]]\n"); 
	}
}# outSameRd() 

sub get_flag {
	my ($ar) = @_; 
	my @keepStr=(); 
	for (my $i=0; $i<@$ar; $i++) {
		push(@keepStr, "$i=$ar->[$i][0]"); 
	}
	@keepStr == 12 or die "@keepStr\n|$keepStr[0]|\n"; 
	my $str = join(",", @keepStr); 
	defined $str2flag{$str} and return($str2flag{$str}); 
	my $ff = &SeqAlnSunhh::mk_flag('keep'=>$str); 
	my @kk = keys %$ff; 
	@kk == 1 or &stopErr("[Err] Bad input string [$str]\n"); 
	$str2flag{$str} = $kk[0]; 
	return($kk[0]); 
}# get_flag() 

