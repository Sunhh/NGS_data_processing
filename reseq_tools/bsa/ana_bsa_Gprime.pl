#!/usr/bin/perl
### Refernece 1 : delta-SNP-index : Takagi, Hiroki, et al. "QTL‐seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations." The Plant Journal 74.1 (2013): 174-183. doi:10.1111/tpj.12105
### Reference 2 : G prime : Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk Segregant Analysis Using Next Generation Sequencing. PLOS Computational Biology 7(11): e1002255. doi: 10.1371/journal.pcbi.1002255
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"windSize:i", # 
	"r_pval:s",   # get_pval_for_Gprime.R 
	"exeRscript:s", # Rscript 
	"cpuN:i",       # 20
	"log_time!", 
); 
$opts{'windSize'} //= 100e3; 
$opts{'exeRscript'} //= 'Rscript'; 
$opts{'cpuN'}       //= 1; 
my $absDir     = &fileSunhh::_dirname( &fileSunhh::_abs_path($0) ); 
my $absProgram = &fileSunhh::_basename( &fileSunhh::_abs_path($0) ); 
my $absCurDir  = &fileSunhh::_abs_path( "./" ); 
unless (defined $opts{'r_pval'}) {
	# get_pval_for_Gprime.R
	$opts{'r_pval'} = "${absDir}/get_pval_for_Gprime.R"; 
}

-e $opts{'r_pval'} or &stopErr("[Err] The Rscript (-r_pval) [$opts{'r_pval'}] is required.\n"); 

my $help_txt = <<HH; 
####################################################################################################
#  perl $0   set1_red2whiteFlesh.bsaTab   -windSize 100000 > set1_red2whiteFlesh.bsaGprime
#
#  -help 
#
#  -windSize     [$opts{'windSize'}] 
#  -r_pval       []
####################################################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

# [Sunhh@bioinfor01 20180521]$ head -5 set1_red2whiteFlesh.bsaTab
# CHROM   POS     n2_RefAD_H      n4_AltAD_H      n1_RefAD_L      n3_AltAD_L
# WK1_R1  73      7       6       16      8
# WK1_R1  78      8       7       16      8
# WK1_R1  170     32      13      23      15
# WK1_R1  174     32      14      23      15

my (@outTbl, @iDp_ar, @Gp_ar); 
my (%forTriCube); # tri-cube smoothing. {'dataType'}{chr} = [ [pos, data_value], [pos, data_value], ... ]; 
my (%order); 
$opts{'log_time'} and &tsmsg("[Msg] Loading input data.\n"); 
while (<>) {
	$order{'num'} ++; 
	chomp; 
	my @ta = &splitL("\t", $_); 
	if ( $ta[0] eq 'CHROM' ) {
		print STDOUT join("\t", @ta, qw/SNPidx_h SNPidx_l SNPidx_d SNPidx_dprime G Gprime nSNPs pvalue negLog10P qvalue negLog10AdjP/)."\n"; 
		next; 
	}
	my ($n1, $n2, $n3, $n4) = @ta[4,2,5,3]; 
	my ($iH, $iL, $iD) = &getHLDIdx( $n1, $n2, $n3, $n4 ); 
	my $G = &getG( $n1, $n2, $n3, $n4 ); 
	push(@outTbl, [ $order{'num'}, \@ta, $iH, $iL, $iD, $G ]); 
	push(@{$forTriCube{'iDprime'}{$ta[0]}}, [ $order{'num'}, $ta[1], $iD ]); 
	push(@{$forTriCube{'Gprime'}{$ta[0]}},  [ $order{'num'}, $ta[1], $G  ]); 
	$order{'chr'}{$ta[0]} //= $order{'num'}; 
}

my %num2val; 
for my $dataType (qw/iDprime Gprime/) {
	my @smoothV; 
	my %idx_smoothV_to_num; 
	if ( $opts{'cpuN'} > 1 ) {
		my $pm  = &LogInforSunhh::get_pm( $opts{'cpuN'} ); 
		my $wd0 = &fileSunhh::new_tmp_dir( 'create'=>1 ); 
		my $sN  = 0; 
		for my $chrID (sort keys %{$forTriCube{$dataType}}) {
			if ( @{$forTriCube{$dataType}{$chrID}} < 5e3 ) {
				# If there are less than 5k sites, I don't need to use mult-threading; 
				$opts{'log_time'} and &tsmsg("[Msg] Getting tri-cube smoothing data for chr[$chrID] of data[$dataType]\n"); 
				my ($v_aref) = &get_TriCube( $forTriCube{$dataType}{$chrID}, $opts{'windSize'} ); 
				my $i_s = scalar(@smoothV) - 1; 
				push(@smoothV, @$v_aref); 
				for (my $i_v=0; $i_v < @$v_aref; $i_v++) {
					$i_s++; 
					$idx_smoothV_to_num{$i_s} = $forTriCube{$dataType}{$chrID}[$i_v][0]; 
				}
				$opts{'log_time'} and &tsmsg("[Msg]   Done for getting tri-cube smoothing data for chr[$chrID] of data[$dataType]\n"); 
				next; 
			}
			$sN ++; 
			my $pid = $pm->start and next; 
			$opts{'log_time'} and &tsmsg("[Msg] Getting tri-cube smoothing data for chr[$chrID] (sN=$sN) of data[$dataType]\n"); 
			my ($v_aref) = &get_TriCube( $forTriCube{$dataType}{$chrID}, $opts{'windSize'} ); 
			my $ofht = &openFH("$wd0/sub_$sN", '>'); 
			for (my $i_v=0; $i_v < @$v_aref; $i_v++) {
				print {$ofht} join("\t", $forTriCube{$dataType}{$chrID}[$i_v][0], @{$v_aref->[$i_v]})."\n"; 
			}
			close($ofht); 
			$opts{'log_time'} and &tsmsg("[Msg]   Done for getting tri-cube smoothing data for chr[$chrID] (sN=$sN) of data[$dataType]\n"); 
			$pm->finish; 
		}
		$pm->wait_all_children; 
		$opts{'log_time'} and &tsmsg("[Msg] Combine sub_files for data[$dataType]\n"); 
		for (my $i=1; $i<=$sN; $i++) {
			my $ifht = &openFH("$wd0/sub_$i", '<'); 
			while (<$ifht>) {
				chomp; 
				my @ta = &splitL("\t", $_); 
				push(@smoothV, [ @ta[1..$#ta] ]); 
				$idx_smoothV_to_num{$#smoothV} = $ta[0]; 
			}
			close($ifht); 
		}
		&fileSunhh::_rmtree($wd0); 
		$opts{'log_time'} and &tsmsg("[Msg]   Done for combining sub_files for data[$dataType]\n"); 
	} else {
		for my $chrID (sort keys %{$forTriCube{$dataType}}) {
			$opts{'log_time'} and &tsmsg("[Msg] Getting tri-cube smoothing data for chr[$chrID] of data[$dataType]\n"); 
			my ($v_aref) = &get_TriCube( $forTriCube{$dataType}{$chrID}, $opts{'windSize'} ); 
			my $i_s = scalar(@smoothV)-1; 
			push(@smoothV, @$v_aref); 
			for (my $i_v=0; $i_v < @$v_aref; $i_v++) {
				$i_s++; 
				$idx_smoothV_to_num{$i_s} = $forTriCube{$dataType}{$chrID}[$i_v][0]; 
			}
			$opts{'log_time'} and  &tsmsg("[Msg]   Done for getting tri-cube smoothing data for chr[$chrID] of data[$dataType]\n"); 
		}
	}
	if ( $dataType eq 'Gprime' ) {
		&add_pval( \@smoothV ); 
	}
	for (my $i=0; $i<@smoothV; $i++) {
		my $t1 = $smoothV[$i]; 
		my $num = $idx_smoothV_to_num{$i}; 
		$num2val{$dataType}{$num} = [@$t1]; 
	}
}


# Output all results; 
$opts{'log_time'} and &tsmsg("[Msg] Output results.\n"); 
for my $t1 (@outTbl) {
	my $num = $t1->[0]; 
	print STDOUT join( "\t", @{$t1->[1]}, @{$t1}[2,3,4], $num2val{'iDprime'}{$num}[0], $t1->[5], @{$num2val{'Gprime'}{$num}} )."\n"; 
}

# Input  : @{$forTriCube{'iDprime'}{$chrID}}, in which each element is [ order_num, pos, data_value ]; 
#            order_num should be unique for each element in array; 
#            This is not a sliding window method; 
# Return : ( [ [$triCube_value, $nsnps], [$triCube_value, $nsnps], ... ] ); 
#          Add tri-cube smoothed value in the same order of the input array; 
sub get_TriCube {
	my ($aref, $wsize) = @_; 
	$wsize //= $opts{'windSize'}; # The window size for each focal SNP site; 
	my $hwsize = ($wsize-1)/2;        # The half size of a window; 
	$hwsize >= 1 or &stopErr("[Err] window size must >= 3 for get_TriCube()\n"); 
	my @srt; # This will be the newly generated sites array with mirrored sites; 
	         # And the 0-element will be replaced by the idx of input array; 
		 # ( [rawIdx, pos, data_value], [rawIdx, pos, data_value], ...); 
	{
		# Replace 0-element with index of input array; 
		my $i = -1; 
		for my $t1 (@$aref) {
			$i++; 
			push(@srt, [$i, $t1->[1], $t1->[2]]); 
		}
		@srt = sort { $a->[1] <=> $b->[1] } @srt; # This will be the newly generated sites array with mirrored sites; 
	}
	my ($i_s, $i_e) = (0, $#srt); # The start/end of real sites in @srt; 
	{
		# Mirror both edges to balance window; 
		my @edge_left  = @{$srt[0]}; 
		my @edge_right = @{$srt[-1]}; 
		MIRROR_RIGHT: 
		for (my $i=$#srt-1; $i>=0; $i--) {
			$srt[$i][1] >= $edge_right[1]-$hwsize or last MIRROR_RIGHT; 
			push( @srt, [ 'NA', 2*$edge_right[1]-$srt[$i][1], $srt[$i][2] ] ); 
		}
		my @toAddLeft; 
		MIRROR_LEFT: 
		for (my $i=1; $i<=$i_e; $i++) {
			$srt[$i][1] <= $edge_left[1]+$hwsize or last MIRROR_LEFT; 
			unshift( @toAddLeft, [ 'NA', 2*$edge_left[1]-$srt[$i][1], $srt[$i][2] ] ); 
		}
		unshift(@srt, @toAddLeft); 
		$i_s += scalar(@toAddLeft); 
		$i_e += scalar(@toAddLeft); 
	}
	my @srt_tricV; # ([$tricube_value, $rawIdx, $nsnps], [$tricube_value, $rawIdx, $nsnps], ...); 
	for (my $i=$i_s; $i<=$i_e; $i++) {
		if ( $srt[$i][2] eq 'NA' ) {
			$srt_tricV[$i] = ['NA', $srt[$i][0], 0]; 
			next; 
		}
		my $curr_pos = $srt[$i][1]; 
		$srt[$i][3] = 1; 
		my $curr_Sw = $srt[$i][3]; 
		my $curr_sum = 0; 
		my ($left_end, $right_end); 
		for (my $le_i=$i-1; $le_i>=0; $le_i--) {
			$srt[$le_i][1] >= $curr_pos-$hwsize or last; 
			$srt[$le_i][2] eq 'NA' and next; 
			my $D = ($curr_pos-$srt[$le_i][1])/$hwsize; 
			$srt[$le_i][3] = (1-$D**3)**3; # Tri-cube; 
			$curr_Sw += $srt[$le_i][3]; 
			$left_end = $le_i; 
		}
		for (my $ri_i=$i+1; $ri_i<@srt; $ri_i++) {
			$srt[$ri_i][1] <= $curr_pos+$hwsize or last; 
			$srt[$ri_i][2] eq 'NA' and next; 
			my $D = ($srt[$ri_i][1]-$curr_pos)/$hwsize; 
			$srt[$ri_i][3] = (1-$D**3)**3; 
			$curr_Sw += $srt[$ri_i][3]; 
			$right_end = $ri_i; 
		}
		$left_end  //= $i; 
		$right_end //= $i; 
		my $nsnp     = 0; # The number of sites used to count tri-cube value; 
		for (my $j=$left_end; $j<=$right_end; $j++) {
			$srt[$j][2] eq 'NA' and next; 
			$nsnp ++; 
			$curr_sum += $srt[$j][2]*$srt[$j][3]/$curr_Sw; 
		}
		push(@srt_tricV, [$curr_sum, $srt[$i][0], $nsnp]); 
	}
	my @back_tricV; 
	my %chk_rawIdx = map { $_=>1 } ( 0 .. $#$aref); 
	for (@srt_tricV) {
		defined $chk_rawIdx{$_->[1]} or &stopErr("[Err] Unknown raw index [$_->[1] : @$_]\n"); 
		$chk_rawIdx{$_->[1]} == 1 or &stopErr("[Err] Repeat raw index [$_->[1] : @$_]\n"); 
		$chk_rawIdx{$_->[1]} ++; 
		$back_tricV[$_->[1]] = [$_->[0], $_->[2]]; 
	}
	for (keys %chk_rawIdx) {
		$chk_rawIdx{$_} == 2 and next; 
		&stopErr("[Err] Bad result for raw index [$_] in raw array [@{$aref->[$_]}]\n"); 
	}
	return(\@back_tricV); 
}# get_TriCube() 


# Retrieve SNP-index for High/Low bulk, and count delta-SNP-index (D=H-L); 
# Return: ($highBulk_SNPindex, $lowBulk_SNPindex, $delta_SNPindex) 
#   Return 'NA' if the value is unavailable. 
# Reference : Takagi, Hiroki, et al. "QTL‐seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations." The Plant Journal 74.1 (2013): 174-183. doi:10.1111/tpj.12105
##########################################
#  .                 Low_bulk  High_bulk
#  A0(ref allele)    n1        n2
#  A1(alt allele)    n3        n4
##########################################
sub getHLDIdx {
	my ($n1, $n2, $n3, $n4) = @_; 
	my $idx_high = ($n2+$n4 > 0) ? $n4/($n2+$n4) : 'NA' ; 
	my $idx_low  = ($n1+$n3 > 0) ? $n3/($n1+$n3) : 'NA' ; 
	my $idx_delt = ( $idx_high ne 'NA' and $idx_low ne 'NA' ) ? $idx_high - $idx_low : 'NA' ; 
	return($idx_high, $idx_low, $idx_delt); 
}# getHLDIdx() 

# Compute G value according to reference_2 ; 
# Return: ($G)
sub getG {
	my ($n1, $n2, $n3, $n4) = @_; 
	my $sum = $n1+$n2+$n3+$n4; 
	$sum > 0 or return('NA');
	my $np1 = ($n1+$n2)*($n1+$n3)/$sum; 
	my $np2 = ($n2+$n1)*($n2+$n4)/$sum; 
	my $np3 = ($n3+$n4)*($n3+$n1)/$sum; 
	my $np4 = ($n4+$n3)*($n4+$n2)/$sum; 
	$np1 > 0 or return('NA'); 
	$np2 > 0 or return('NA'); 
	$np3 > 0 or return('NA'); 
	$np4 > 0 or return('NA'); 
	my $v1 = ( $n1 > 0 ) ? $n1 * log($n1/$np1) : 0; 
	my $v2 = ( $n2 > 0 ) ? $n2 * log($n2/$np2) : 0; 
	my $v3 = ( $n3 > 0 ) ? $n3 * log($n3/$np3) : 0; 
	my $v4 = ( $n4 > 0 ) ? $n4 * log($n4/$np4) : 0; 
	my $g = 2 * ($v1+$v2+$v3+$v4); 
	return($g); 
}# getG() 

# Count p-value for iDp and Gp; 
# Input  : ( $tricubeV_aref ) 
#            $tricubeV_aref is output of &get_TriCube(); 
# Return : ( [ [pvalue, negLog10P, qvalue_fdr, negLog10AdjP], [pvalue, negLog10P, qvalue_fdr, negLog10AdjP], ... ] )
# Return : undef(); Add [pvalue, negLog10P, qvalue_fdr, negLog10AdjP] to each element of input $tricubeV_aref; 
sub add_pval {
	my ($aref) = @_; 
	my $wd = &fileSunhh::new_tmp_dir( 'create'=>1 ); 
	&tsmsg("[Msg] Working for p-value in dir [$wd]\n"); 
	my $ofh = &openFH("$wd/inF", '>'); 
	for (@$aref) {
		print {$ofh} "$_->[0]\n"; 
	}
	close($ofh); 
	&exeCmd_1cmd("$opts{'exeRscript'} $opts{'r_pval'} $wd/inF >/dev/null"); 
	my $ifh = &openFH("$wd/inF.pval", '<'); 
	my @back; 
	while (<$ifh>) {
		chomp; 
		my @ta = &splitL("\t", $_); 
		$. == 1 and next; 
		# push(@back, [@ta[0,1,2,3]]); 
		push(@back, [@ta]); 
	}
	close($ifh); 
	&fileSunhh::_rmtree($wd); 
	scalar(@back) == scalar(@{$aref}) or &stopErr("[Err] Input and output number different!\n"); 
	for (my $i=0; $i<@back; $i++) {
		push(@{$aref->[$i]}, @{$back[$i]}); 
	}
	return(\@back); 
}# add_pval() 



