#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"libInfoFn:s", 
	"eleInfoFn:s", 
	"rdCntFn:s", 
	"skipEleID:s", 
	"cntRPKM:s", 
	"help!", 
); 

my %gg; 
&set_glob(); 

$gg{'libInfoFn'} ne '' and &load_libInfo(); 
$gg{'skipEleID'} ne '' and &load_skipEleID(); 

&load_eleInfo(); 

&load_rdCnt(); 

&count_tpm(); 

for (@{$gg{'tpmTbl'}}) {
	print STDOUT join("\t", @$_)."\n"; 
}

if (defined $gg{'cntRPKM'}) {
	my $ofh = &openFH("$gg{'cntRPKM'}", '>'); 
	for (@{$gg{'rpkmTbl'}}) {
		print {$ofh} join("\t", @$_)."\n"; 
	}
	close($ofh); 
}

######################################################################################
#  Sub-routines
######################################################################################
sub set_glob {
	# $gg{'libRdLen'} = 100; 
	$gg{'help_txt'} = <<"H1"; 
################################################################################
# perl $0 -rdCntFn rdCnt_SenseStrand   -eleInfoFn gen_len_file
#
# -help
#
# -rdCntFn           [filename] Required. Format : EleID \\t lib1_Cnt \\t lib2_Cnt
# -eleInfoFn         [filename] Required. 
#                        'trans_len' is used as transcript legnth to compute TPM; 
#                        Format : 
#                          eleID \\t value_type \\t value
#                          ele_1 \\t trans_len  \\t 1035
#                          ele_2 \\t trans_len  \\t 999
#
# -skipEleID         [filename] A file listing the element IDs not to be included in any analysis. 
#                        For example : __no_feature / __ambiguous that output from htseq-count; 
#
# -libInfoFn         [filename] Used to add/replace information. 
#                        'lib_tNSum_woRdLen' will replace the result computed from -eleInfoFn; 
#                        Format : 
#                          libID \\t value_type  \\t value
#                          lib1  \\t lib_tNSum_woRdLen \\t 88888823
#
#
# -cntRPKM           [filename] A file storing RPKM values. 
#
################################################################################
H1
#                          lib1  \\t lib_rd_len  \\t 100    # This is not used currently, as I assume all the reads in one library has the same read length. 
# -libRdLen          [$gg{'libRdLen'}] read length in each library by default. 
#                        This is not used currently, as I assume all the reads in one library has the same read length. 
#                        If the read length in libraries are different, please assign them in -libInfoFn
#                        For paired-end sequencing, this value should be also single-read-length. 
	$opts{'help'} and &LogInforSunhh::usage($gg{'help_txt'}); 
	defined $opts{'rdCntFn'} or &stopErr("[Err] -rdCntFn is required.\n"); 
	defined $opts{'eleInfoFn'} or &stopErr("[Err] -eleInfoFn is required.\n"); 
	# for my $k1 (qw/rdCntFn eleInfoFn libRdLen/) {
	for my $k1 (qw/rdCntFn eleInfoFn cntRPKM/) {
		defined $opts{$k1} and $gg{$k1} = $opts{$k1}; 
	}
	$gg{'libInfoFn'} = ( defined $opts{'libInfoFn'} ) ? $opts{'libInfoFn'} : ''; 
	$gg{'skipEleID'} = ( defined $opts{'skipEleID'} ) ? $opts{'skipEleID'} : ''; 
	return; 
}# set_glob() 

sub load_libInfo {
	my $fh = &openFH($gg{'libInfoFn'}, '<'); 
	while (<$fh>) {
		m!^\s*(#|$)! and next; 
		chomp; 
		my @ta = &splitL("\t", $_); 
		$gg{'libInfo'}{$ta[0]}{$ta[1]} = $ta[2]; 
	}
	close($fh); 
	return; 
}# load_libInfo() 

sub load_skipEleID {
	my $fh = &openFH($gg{'skipEleID'}, '<'); 
	while (<$fh>) {
		chomp; 
		my @ta = &splitL("\t", $_); 
		$gg{'skipEleID_hash'}{$ta[0]} = 1; 
	}
	close($fh); 
	return; 
}# load_skipEleID

sub load_eleInfo {
	my $fh = &openFH($gg{'eleInfoFn'}, '<'); 
	while (<$fh>) {
		m!^\s*(#|$)! and next; 
		my @ta = &splitL("\t", $_); 
		$gg{'eleInfo'}{$ta[0]}{$ta[1]} = $ta[2]; 
	}
	close($fh); 
}# load_eleInfo() 

sub load_rdCnt {
	my $fh = &openFH($gg{'rdCntFn'}, '<'); 
	my $h_txt = <$fh>; 
	chomp($h_txt); 
	$gg{'cntTbl'}[0] = [ &splitL("\t", $h_txt) ]; 
	for (my $i=1; $i<@{$gg{'cntTbl'}[0]}; $i++) {
		$gg{'lib2idx'}{$gg{'cntTbl'}[0][$i]} = $i; 
		$gg{'idx2lib'}{$i} = $gg{'cntTbl'}[0][$i]; 
	}
	while (<$fh>) {
		chomp; 
		my @ta = &splitL("\t", $_); 
		defined $gg{'skipEleID_hash'}{$ta[0]} and next; 
		push(@{$gg{'cntTbl'}}, [@ta]); 
	}
	close($fh); 
}# load_rdCnt() 

sub count_tpm {
	# Firstly, I need to check and complete any missing information for computing. 
	for my $libID (@{$gg{'cntTbl'}[0]}[ 1 .. $#{$gg{'cntTbl'}[0]} ]) {
		# $gg{'libInfo'}{$libID}{'lib_rd_len'} //= $gg{'libRdLen'}; 
		$gg{'libInfo'}{$libID}{'lib_tNSum_woRdLen'} //= -1; # -1 indicates I need to count. 
		defined $gg{'cntRPKM'} and $gg{'libInfo'}{$libID}{'lib_sumMapRd'} //= -1; # 
	}
	for my $eleID (map { $_->[0] } @{$gg{'cntTbl'}}[1 .. $#{$gg{'cntTbl'}}] ) {
		defined $gg{'eleInfo'}{$eleID}{'trans_len'} or &stopErr("[Err] Failed to find trans_len for eleID [$eleID]\n"); 
	}
	# Then I need to count each gene's TPM without nomalization; 
	my @lib_tN_est1; # Not divided by read_length; And it is 1e6 times of lib_tNSum_woRdLen; 
	my @lib_sumMapRd; # This is used for RPKM computation; 
	for (my $i=1; $i<@{$gg{'cntTbl'}}; $i++) {
		my $ar1 = $gg{'cntTbl'}[$i]; 
		my $eleID = $ar1->[0]; 
		# trans_len : $gg{'eleInfo'}{$eleID}{'trans_len'}; 
		for (my $j=1; $j<@$ar1; $j++) {
			my $libID = $gg{'cntTbl'}[0][$j]; 
			# read_len  : $gg{'libInfo'}{$libID}{'lib_rd_len'}; 
			# cuntRd    : $ar1->[$j]
			my $ele_tN_est1 = $ar1->[$j] * 1e6 / $gg{'eleInfo'}{$eleID}{'trans_len'}; # (rdNum * 1e6) / gene_length ; Not dvd by read_length; 
			# I want to re-caculate these values in order to reduce approximate computation. 
			$lib_tN_est1[$j] += $ele_tN_est1; 
			defined $gg{'cntRPKM'} and $lib_sumMapRd[$j] += $ar1->[$j]; # For RPKM 
		}
	}
	# Count lib_tNSum_woRdLen for each library; 
	for (my $j=1; $j<@{$gg{'cntTbl'}[0]}; $j++) {
		my $libID = $gg{'cntTbl'}[0][$j]; 
		if ( $gg{'libInfo'}{$libID}{'lib_tNSum_woRdLen'} < 0 ) {
			$gg{'libInfo'}{$libID}{'lib_tNSum_woRdLen'} = $lib_tN_est1[$j]/1e6; 
		} 
		if ( $gg{'libInfo'}{$libID}{'lib_tNSum_woRdLen'} == 0 ) {
			&tsmsg("[Wrn] Bad 'lib_tNSum_woRdLen' for lib [$libID] with value [$gg{'libInfo'}{$libID}{'lib_tNSum_woRdLen'}] changed to 1\n"); 
			$gg{'libInfo'}{$libID}{'lib_tNSum_woRdLen'} = 1; 
		}
		if ( defined $gg{'cntRPKM'} ) {
			# For RPKM 
			if ( $gg{'libInfo'}{$libID}{'lib_sumMapRd'} < 0 ) {
				$gg{'libInfo'}{$libID}{'lib_sumMapRd'} = $lib_sumMapRd[$j]; 
			} elsif ( $gg{'libInfo'}{$libID}{'lib_sumMapRd'} == 0 ) {
				&tsmsg("[Wrn] Bad 'lib_sumMapRd' for lib [$libID] with value [$gg{'libInfo'}{$libID}{'lib_sumMapRd'}] changed to 1\n"); 
				$gg{'libInfo'}{$libID}{'lib_sumMapRd'} = 1; 
			}
		}
	}
	# Normalize read counts to TPM; 
	$gg{'tpmTbl'}[0] = [ @{$gg{'cntTbl'}[0]} ]; 
	defined $gg{'cntRPKM'} and $gg{'rpkmTbl'}[0] = [ @{$gg{'cntTbl'}[0]} ]; 
	for (my $i=1; $i<@{$gg{'cntTbl'}}; $i++) {
		my $ar1 = $gg{'cntTbl'}[$i]; 
		my $eleID = $ar1->[0]; 
		$gg{'tpmTbl'}[$i][0] = $eleID; 
		defined $gg{'cntRPKM'} and $gg{'rpkmTbl'}[$i][0] = $eleID; 
		# trans_len : $gg{'eleInfo'}{$eleID}{'trans_len'}; 
		for (my $j=1; $j<@$ar1; $j++) {
			my $libID = $gg{'cntTbl'}[0][$j]; 
			# read_len  : $gg{'libInfo'}{$libID}{'lib_rd_len'}; 
			# cuntRd    : $ar1->[$j]
			$gg{'tpmTbl'}[$i][$j] = $ar1->[$j] * 1e6 / ( $gg{'eleInfo'}{$eleID}{'trans_len'} * $gg{'libInfo'}{$libID}{'lib_tNSum_woRdLen'} ); 
			defined $gg{'cntRPKM'} and $gg{'rpkmTbl'}[$i][$j] = $ar1->[$j] * 1e9 / ( $gg{'eleInfo'}{$eleID}{'trans_len'} * $gg{'libInfo'}{$libID}{'lib_sumMapRd'} ); 
		}
	}
	return; 
}# count_tpm() 




