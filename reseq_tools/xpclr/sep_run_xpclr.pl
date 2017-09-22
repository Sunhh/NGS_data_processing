#!/usr/bin/perl -w
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

!@ARGV and die "perl $0 in_prefix xpclr_para sep_siteN out_xpclr\n"; 

my $in_pref = shift;
my $xpclr_para = shift; 
my $siteN = shift; 
my $ofile = shift; 

my %glob; 
$glob{'fh_out'} = &openFH($ofile,'>'); 

my %xpclr_hash = %{ &parse_para($xpclr_para) }; 

my @geno_1 = @{ &load_geno("${in_pref}_g1.geno") }; 
my @geno_2 = @{ &load_geno("${in_pref}_g2.geno") }; 
my @snp    = @{ &load_snp("${in_pref}.snp") }; 

my %sep; 
my %cnt; 
my $tmp_dir = &fileSunhh::new_tmp_dir(); 
mkdir($tmp_dir) or &stopErr("[Err] Failed to create dir [$tmp_dir]\n"); 

# Look for the effective start P
$glob{'ladder'} = [200 , 2000, 20000]; 
for my $n (@{$glob{'ladder'}}) {
	my $tn = $n; 
	$tn > $#snp and $tn = $#snp; 
	if ( $tn == -1 ) {
		# There is no input at all. 
		exit; 
	}
	&fileSunhh::write2file( "${tmp_dir}/test_g1.geno", join("\n", @geno_1[0 .. ($tn-1)])."\n", '>' ); 
	&fileSunhh::write2file( "${tmp_dir}/test_g2.geno", join("\n", @geno_2[0 .. ($tn-1)])."\n", '>' ); 
	&fileSunhh::write2file( "${tmp_dir}/test.snp", join("\n", map { $_->[1] } @snp[0 .. ($tn-1)])."\n", '>' ); 
	&exeCmd_1cmd("XPCLR -xpclr   ${tmp_dir}/test_g1.geno ${tmp_dir}/test_g2.geno ${tmp_dir}/test.snp ${tmp_dir}/test.out   $xpclr_para"); 
	open F,'<',"${tmp_dir}/test.out.xpclr.txt" or die ;
	while (<F>) {
		chomp; 
		my @ta = split(/\s+/, $_); 
		$ta[3] =~ s!\.0+$!! or die "Bad [$ta[3]]\n"; 
		$ta[3] =~ m!^\d+$! or die "Bad $_\n"; 
		$glob{'first_snp_pos'} = $ta[3]; 
		&tsmsg("[Msg] Found first_snp_pos [$glob{'first_snp_pos'}]\n"); 
		last; 
	}
	close F; 
	defined $glob{'first_snp_pos'} and last; 
}
defined $glob{'first_snp_pos'} or &stopErr("[Err] Failed to find first position.\n"); 

for (my $i=0; $i<@snp; $i++) { 
	$snp[$i][0] >= $glob{'first_snp_pos'} or next; 
	$snp[$i][0] == $glob{'first_snp_pos'} or &stopErr("[Err] Bad first_snp_pos [$glob{'first_snp_pos'}] found.\n"); 
	@snp    = @snp[ $i .. $#snp ];
	@geno_1 = @geno_1[ $i .. $#geno_1 ]; 
	@geno_2 = @geno_2[ $i .. $#geno_2 ]; 
	unlink("$tmp_dir/test_g1.geno"); 
	unlink("$tmp_dir/test_g2.geno"); 
	unlink("$tmp_dir/test.snp"); 
	unlink("$tmp_dir/test.out.xpclr.txt");
	unlink("$tmp_dir/test.out.xpclr.log");
	last; 
}


# Prepare separate input. 
#   It seems 400k snp sites is the maximum number of input. 
#   There are very small differences for each grid's xpclr_score with gridSize 20000, which are smaller than 0.0001; 
#   But when using the real data to compare, the xpclr_score with gridSize differs according to the site number used in whole dataset and separated dataset. 
#       If sep_gridSize is 10% of whole_gridSize, xpclr_score differs around 0.1 ; 
#   When using different datasets, the estimated 'w' is different, and 'w' is used in logPRatio (xpclr_score) calculation. 
$cnt{'curr_ID'}   = 1; 
$cnt{'curr_si'}   = 0; 
$sep{$cnt{'curr_ID'}}{'pos_start'} = $snp[$cnt{'curr_si'}][0]; 
&renew_cnt_hash(); 
for (my $i=1; $i<@snp; $i++) {
	if ( $cnt{'curr_siteN'} >= $siteN ) {
		for my $fh ( qw/curr_fh_g1 curr_fh_g2 curr_fh_snp/ ) {
			close ($fh); 
		}
		$i == $#snp and last; 
		$cnt{'curr_ID'} ++; 
		&renew_cnt_hash(); 
		my $prev_position1 = $snp[$cnt{'prev_ei'}][0] - $xpclr_hash{'gridSize'} * 3; 
		my $prev_gm1       = $snp[$cnt{'prev_ei'}][2] - $xpclr_hash{'gWin'} * 3; 
		for (my $j=$cnt{'prev_ei'}; $snp[$j][0] > $prev_position1 or $snp[$j][2] > $prev_gm1; $j--) {
			$j > $cnt{'prev_si'} or &stopErr("[Err] Failed to find a good i at [$j]\n"); 
			$i = $j-1; 
		}
		# Now the position $snp[$i][0] is smaller than ($snp[ $cnt{'prev_ei'} ][0] - $xpclr_hash{'gridSize'}).  
		$cnt{'curr_si'} = $i; 
		$sep{$cnt{'curr_ID'}}{'pos_start'} = $snp[$cnt{'curr_si'}][0]; 
	}
	$cnt{'curr_ei'} = $i; 
	$sep{$cnt{'curr_ID'}}{'pos_end'} = $snp[ $cnt{'curr_ei'} ][0]; 
	$cnt{'curr_siteN'} ++; 
	print {$cnt{'curr_fh_g1'}}  $geno_1[$i]."\n"; 
	print {$cnt{'curr_fh_g2'}}  $geno_2[$i]."\n"; 
	print {$cnt{'curr_fh_snp'}} $snp[$i][1]."\n"; 
}

my %out; 
for my $id (sort {$a<=>$b} keys %sep) {
	my $pref = $sep{$id}{'pref'}; 
	&exeCmd_1cmd("XPCLR -xpclr   ${pref}_g1.geno ${pref}_g2.geno ${pref}.snp ${pref}   $xpclr_para"); 
	open F,'<',"${pref}.xpclr.txt" or &stopErr("[Err] Failed to find ${pref}.xpclr.txt\n"); 
	while (<F>) {
		chomp; 
		my @ta=split(/ /, $_); 
		$#ta == 6 or &stopErr("[Err] Bad line: $_\n"); 
		my $tk = "$ta[0] $ta[1]"; 
		if ( defined $out{$tk} ) {
			$out{$tk}[3] eq $ta[3] or &stopErr("[Err] Diff lines: \nPrev @{$out{$tk}}\nCurr @ta\n"); 
			$out{$tk}[2] < $ta[2] and $out{$tk} = [@ta]; 
		} else {
			$out{$tk} = [@ta]; 
		}
	}
	close(F); 
}
for my $tk (sort { $out{$a}[1] <=> $out{$b}[1] } keys %out) {
	print {$glob{'fh_out'}} join(' ', @{$out{$tk}})."\n"; 
}
close($glob{'fh_out'}); 

&fileSunhh::_rmtree($tmp_dir); 

&tsmsg("[Rec] All done [$0]\n"); 

sub renew_cnt_hash {
	$cnt{'curr_pref'} = "$tmp_dir/$cnt{'curr_ID'}"; 
	&tsmsg("[Msg] Processing $cnt{'curr_pref'}\n"); 
	# $cnt{'curr_s_posi'} = 0; 
	$cnt{'curr_fh_g1'} = &openFH("$cnt{'curr_pref'}_g1.geno", '>'); 
	$cnt{'curr_fh_g2'} = &openFH("$cnt{'curr_pref'}_g2.geno", '>'); 
	$cnt{'curr_fh_snp'}= &openFH("$cnt{'curr_pref'}.snp",     '>'); 
	$cnt{'curr_si'} //= 0; 
	$cnt{'curr_ei'} //= 0; 
	$cnt{'prev_si'} = $cnt{'curr_si'}; 
	$cnt{'prev_ei'} = $cnt{'curr_ei'}; 
	print {$cnt{'curr_fh_g1'}}  $geno_1[0]."\n"; 
	print {$cnt{'curr_fh_g2'}}  $geno_2[0]."\n"; 
	print {$cnt{'curr_fh_snp'}} $snp[0][1]."\n"; 
	$sep{ $cnt{'curr_ID'} }{'ID'}   = $cnt{'curr_ID'}; 
	$sep{ $cnt{'curr_ID'} }{'pref'} = $cnt{'curr_pref'}; 
	$cnt{'curr_siteN'}= 1; 
}
sub parse_para {
	my $p = shift; 
	$p =~ m/^\s*\-(w\d)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\-(p\d+)\s+(\S+)\s*$/i or &stopErr("[Err] Bad para [$p]\n"); 
	#             -w1 0.0005 100 100 1 -p0 0.7
	my %back; 
	@back{qw/w1 gWin snpWin gridSize chrN p0 corrLevel/} = ($1, $2,   $3, $4, $5,$6,$7); 
	return(\%back); 
}


sub load_snp {
	my $fn = shift; 
	my $fh = &openFH($fn,'<'); 
	my @back; 
	while (<$fh>) {
		chomp; 
		my @ta=split(/\t/, $_);
		push(@back, [$ta[3], $_, $ta[2]]); 
	}
	close($fh); 
	return(\@back); 
}

sub load_geno {
	my $fn = shift; 
	my $fh = &openFH($fn,'<'); 
	my @back; 
	while (<$fh>) {
		chomp; 
		push(@back, $_); 
	}
	close($fh); 
	return(\@back); 
}
