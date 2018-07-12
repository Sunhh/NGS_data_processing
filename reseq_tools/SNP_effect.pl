#!/usr/bin/perl -w 
# Updated on 2018/07/10 Try to accept 
# Updated on 2018/05/24
# Method 1 : 20120110 I don't care about the effect of SNPs combination, so I only compute the changes caused by the independent SNP site; 
#            Classes of SNPs : Only two alleles considered. 
#              1. Intergenic : 
#                 1.1 In the upstream XX-kb of gene; 
#                 1.2 In the downstream XX-kb of gene; 
#                 1.3 Other intergenic region; 
#              2. In intron : 
#                 2.1 In the two bases of a intron's end; (Which may change the splicing); 
#                 2.2 Other intron region; 
#              3. Coding region : 
#                 3.1 In the first codon, and change a start codon to a non-start codon; 
#                 3.2 In a stop codon, and change it to a non-stop codon; 
#                 3.3 Change the AA codon. 
#                 3.4 Change a AA codon to a stop codon; 
#                 3.5 Doesn't change the AA coding; 
#                 3.6 Change a stop codon to other stop codon; 
#                 3.7 Change an initial start codon to other start codon; 
#                 3.8 There is a N or someother problem in reference or sample codon; 
#                 3.9 InDel related coding region; (20180503)
# 方法一: 20120110 不考虑邻近SNP之间的互相影响，独立检查SNP所带来的变化;
#         SNP影响分类: 
#            1. 基因间区 : 1.1 基因上游500bp 1.2 基因下游500bp 1.3 其余基因间区; 
#            2. 基因内含子区: 2.1 intron/exon边界(即exon边界碱基外延2bp) 2.2 其它内含子区域; 
#            3. Coding区: 3.1 起始密码子变为非起始密码子  3.2 终止密码子变化为继续编码  
#                         3.3 氨基酸编码变化-转义  3.4 氨基酸编码变化-终止   3.5 同义突变; 
#                         3.6 终止密码子替换为其它终止密码子; 3.7 起始密码子变为其它起始密码(AA变化); 
# edit20120118:           3.8 位于编码区, 但氨基酸编码中存在'N', 因而无法判断; 
# ClV6 genome intron length stat: SUM             MEAN                    MEDIAN  MIN     MAX     Count   NoNull
#                                 38918967        464.338157392383        199     11      9976    83816   83816
# 2013-08-19 开始遇到太多Indel、杂合等问题，这些单独处理太麻烦，而且有时基因位置信息也是有用的。

use strict; 
use warnings; 
use LogInforSunhh; 
use gffSunhh;
use fastaSunhh;
use fileSunhh;
use mathSunhh;
use SNP_tbl; 
use Getopt::Long; 
my %opts; 

GetOptions(\%opts,
	"cdsGff:s", # Required; 
	"cdsFas:s", "chrFas:s", # At least one is required. 
	"in_snp:s", # Input SNP table; 
	"outEff:s", # 
	"upstreamLen:i", 
	"downstreamLen:i", 
	"asSnpCol!", 
	"help!", 
); 

$opts{'upstreamLen'} //= 500; 
$opts{'downstreamLen'} //= 500; 

my $help_txt = <<HH; 
######################################################################
# perl $0 -cdsGff cds_loc.gff -chrFas chr_seq.fasta -in_snp in_SNP.cols
# 
# -cdsFas [filename] This is required if -chrFas is not provided; 
#
# -cdsGff [Filename] There must be CDS feature in it. 
#
# -upstreamLen    [$opts{'upstreamLen'}]
# -downstreamLen  [$opts{'downstreamLen'}]
#
# -outEff         [Filename] STDOUT instead 
#
# -help   [Boolean]
######################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
defined $opts{'in_snp'} or &LogInforSunhh::usage($help_txt); 
defined $opts{'cdsGff'} or &LogInforSunhh::usage($help_txt); 
( defined $opts{'chrFas'} or defined $opts{'cdsFas'} ) or &LogInforSunhh::usage($help_txt); 

my $gff_obj = gffSunhh->new(); 
my $fas_obj = fastaSunhh->new(); 
my $mat_obj = mathSunhh->new(); 

my $gffFh = &openFH($opts{'cdsGff'}, '<'); 
my $snpFh = &openFH($opts{'in_snp'}, '<'); 
my $outFh = \*STDOUT; 
defined $opts{'outEff'} and $outFh = &openFH($opts{'outEff'}, '>'); 
my $wrk_dir = &fileSunhh::new_tmp_dir('create'=>1); 

my %glob; 
$glob{'file_cds2Scf_agp'} = "$wrk_dir/cds2Scf.agp"; 
$glob{'fh_cds2Scf_agp'} = &openFH($glob{'file_cds2Scf_agp'}, '>'); 
$glob{'file_intron2Scf_agp'} = "$wrk_dir/intron2Scf.agp"; 
$glob{'fh_intron2Scf_agp'} = &openFH($glob{'file_intron2Scf_agp'}, '>'); 
$glob{'file_up2Scf_agp'} = "$wrk_dir/upstream2Scf.agp"; 
$glob{'fh_up2Scf_agp'} = &openFH($glob{'file_up2Scf_agp'}, '>'); 
$glob{'file_down2Scf_agp'} = "$wrk_dir/downstream2Scf.agp"; 
$glob{'fh_down2Scf_agp'} = &openFH($glob{'file_down2Scf_agp'}, '>'); 
$glob{'d2b'} = { &SNP_tbl::get_diploid_d2b() }; 
for my $t0 (qw/A T G C/) {
	$glob{'d2b'}{$t0} = [$t0]; 
}
for (keys %opts) {
	$glob{$_} //= $opts{$_}; 
}


##### Loading gff file and chr fasta seq file; 
$glob{'gff_top_hier'} = { 'mrna'=>1, 'match'=>2, 'protein_match'=>3, 'expressed_sequence_match'=>4 }; 
my (%in_gff, %in_seq); 
&tsmsg("[Msg] Loading gff [$opts{'cdsGff'}]\n"); 
&load_gff($gffFh, \%in_gff, \%in_seq); 
if ( defined $opts{'chrFas'} ) {
	&tsmsg("[Msg] Loading chr fasta [$opts{'chrFas'}]\n"); 
	%in_seq = %{ $fas_obj->save_seq_to_hash('faFile'=>$opts{'chrFas'}) }; 
}
for (keys %in_seq) {
	$in_seq{$_}{'seq'} =~ s!\s!!g; 
	$in_seq{$_}{'seq'} = uc($in_seq{$_}{'seq'}); 
	$in_seq{$_}{'len'} = length($in_seq{$_}{'seq'}); 
}
##### Loading fasta files and generate CDS sequences; 
my (%in_cds, %in_intron); # frame is one of [1,2,3]
if ( defined $opts{'cdsFas'} ) {
	&tsmsg("[Msg] Loading CDS fasta [$opts{'cdsFas'}]\n"); 
	%in_cds = %{ $fas_obj->save_seq_to_hash('faFile'=>$opts{'cdsFas'}) }; 
	for my $k (keys %in_cds) {
		$in_cds{$k}{'seq'} =~ s!\s!!g; 
		$in_cds{$k}{'len'} = length($in_cds{$k}{'seq'}); 
		if ($in_cds{$k}{'definition'} =~ m!\[frame=(\d+)\]!i) {
			$in_cds{$k}{'frame'} = $1; 
		} else {
			$in_cds{$k}{'frame'} = 1; 
		}
	}
}
my (%agp_cds2Scf, %agp_intron2Scf, %agp_up2Scf, %agp_down2Scf); 
&tsmsg("[Msg] Retriving CDS from gff\n"); 
&retrive_cds(); 
for (keys %in_cds) {
	$in_cds{$_}{'seq'} =~ s!\s!!g; 
	$in_cds{$_}{'seq'} = uc($in_cds{$_}{'seq'}); 
	$in_cds{$_}{'len'} = length($in_cds{$_}{'seq'}); 
}

##### All information is gotten, so I can start to process the snp table. 
&tsmsg("[Msg] Processing SNPs\n"); 
SNP_LINE:
while (my $l = &wantLineC($snpFh)) {
	$l =~ m!^\s*$! and next; 
	my @ta = &splitL("\t", $l); 
	my ($chr_id, $chr_pos) = @ta[0,1]; 
	if ( $chr_id =~ m!^#?(chr|chrom|chomosome)$!i ) {
		print {$outFh} join("\t", "SNP_Effect", $l)."\n"; 
		next SNP_LINE; 
	}
	
	# Read in each allele from left to right; 
	my $ref_allele; 
	my @alt_allele; 
	my @effects; 
	my %used_allele; 
	my @pop_allele; 
	
	# Process each genotype to store alleles; 
	for my $tb (@ta[2..$#ta]) {
		$tb eq './.' and next; 
		$tb =~ m!^n*$!i and next; 
		$tb = uc($tb); 
		my @tc; # This stores all alleles in current genotype: A/T/G/C/Del/Ins/InDel
		if ($tb =~ m!^([ATGC\*N]+)/([ATGC\*N]+)$!) {
			@tc = ($1, $2); 
			push(@pop_allele, @tc); 
			if (length($1) > 1 or length($2) > 1) {
				@tc = ('InDel', 'InDel'); 
			}
		} elsif (defined $glob{'d2b'}{$tb}) {
			@tc = (@{$glob{'d2b'}{$tb}}); 
			push(@pop_allele, @tc); 
		} elsif ($tb eq '*' or $tb eq '-') {
			push(@pop_allele, $tb); 
			@tc = ('Del'); 
		} elsif ($tb =~ m!\+!) {
			# In fact, m!^[^+]++! means a heterozygous insertion, but I don't want it too complex. 
			push(@pop_allele, $tb); 
			@tc = ('Ins'); 
		} elsif ( $tb =~ m!^[ATGCN]{2,}$! ) {
			if ( $opts{'asSnpCol'} ) {
				if ( $tb =~ m!^([ATGC])([ATGC])$! ) {
					@tc = ($1, $2); 
					push(@pop_allele, @tc); 
				} else {
					&tsmsg("[Wrn] Skip unknown genotype [$tb]\n"); 
					next; 
				}
			} else {
				push(@pop_allele, $tb); 
				@tc = ('InDel', 'InDel'); 
			}
		} else {
			&tsmsg("[Wrn] Skip unknown genotype [$tb]\n"); 
			next; 
		}
		# Store reference and alternative alleles; 
		@tc > 0 or next; 
		unless ( defined $ref_allele ) {
			for my $tc0 (@tc) {
				$tc0 =~ m!^[ATGC]$! or next; 
				$ref_allele = $tc0; 
				$used_allele{$ref_allele} = 1; 
			}
		}
		for my $td (@tc) {
			defined $used_allele{$td} and next; 
			push(@alt_allele, $td); 
			$used_allele{$td} = 1; 
		}
	}# for my $tb (@ta[2..$#ta]) : 
	unless ( defined $ref_allele ) {
		for my $tp0 (@pop_allele) {
			defined $used_allele{$tp0} and next; 
			$tp0 =~ m!^[ATGC]$! and do { $ref_allele = $tp0; last; }; 
		}
		$ref_allele //= 'N'; 
		$used_allele{$ref_allele} = 1; 
	}
	
	# Check the effects of alleles; 
	die "|@alt_allele|\n"; 
	for (my $i=0; $i<@alt_allele; $i++) {
		my @t_eff = &snp_eff($ref_allele, $alt_allele[$i], $chr_id, $chr_pos, \%agp_cds2Scf, \%agp_intron2Scf, \%agp_up2Scf, \%agp_down2Scf); 
		push(@effects, join(",", map { join(":", @$_) } @t_eff)); 
	}
	@effects == 0 and push(@effects, "Skipped"); 
	
	print {$outFh} join("\t", $effects[0], $l)."\n"; 
}# End while ()
close($snpFh); 
&fileSunhh::_rmtree($wrk_dir); 

######################################################################
# Subroutines : 
######################################################################

sub snp_eff {
	my ($rb, $ab, $chr_id, $chr_pos, $agp_c2s_href, $agp_i2s_href, $agp_u2s_href, $agp_d2s_href) = @_; 
	my ($rb_ori, $ab_ori) = ($rb, $ab); 
	my @back_eff; 
	# For cds2scf : Check if position is in CDS; 
	my @cds_eff; 
	my @got_c2s = $mat_obj->switch_position(
		'qry2ref' => $agp_c2s_href, 
		'qryID'   => $chr_id, 
		'qryPos'  => $chr_pos, 
		'qryStr'  => '+'
	); 
	if (defined $got_c2s[0][0]) {
		# This means this position sits in at least one gene's CDS; 
		for my $cp (@got_c2s) {
			my ($cds_id, $cds_pos, $cds_str) = @$cp; 
			# Get the AA position; 
			$in_cds{$cds_id}{'frame'} //= 1; 
			my $aa_pos   = int(($cds_pos+$in_cds{$cds_id}{'frame'}-2)/3)+1; 
			my $aa_frame = ($cds_pos+$in_cds{$cds_id}{'frame'}-2)%3+1; # The frame is one of [1,2,3]; 
			my $codon_S  = $cds_pos-$aa_frame+1; 
			my $codon_E  = $cds_pos-$aa_frame+3; 
			
			$rb = $rb_ori; 
			$ab = $ab_ori; 
			if ( $rb =~ m!^(InDel|Del|Ins)$!i or $ab =~ m!^(InDel|Del|Ins)$!i ) {
				push(@cds_eff, [$cds_id, '3.9', "${aa_pos}_InDel"]); 
				next; 
			}
			if ($cds_str eq '-') {
				# We need to get the reverse complement; 
				($rb) = &rev_comp($rb); 
				($ab) = &rev_comp($ab); 
			}
			if ($codon_S <= 0 or $codon_E > $in_cds{$cds_id}{'len'}) {
				push(@cds_eff, [$cds_id, '3.8', "X${aa_pos}X"]); 
				next; 
			}
			my $refb3 = substr($in_cds{$cds_id}{'seq'}, $codon_S-1, 3); 
			my $altb3 = $refb3; 
			substr($refb3, $aa_frame-1, 1) = $rb; # The real reference bbb; 
			substr($altb3, $aa_frame-1, 1) = $ab; # The real alternative bbb; 
			my ($refAA, $refStart) = &fastaSunhh::bbb2aa( $refb3, 1 ); 
			my ($altAA, $altStart) = &fastaSunhh::bbb2aa( $altb3, 1 ); 
			if ($refAA eq 'X' or $altAA eq 'X') {
				push(@cds_eff, [$cds_id, '3.8', "X${aa_pos}X"]); 
				next; 
			}
			if ($aa_pos == 1 and $refStart == 1) {
				# Changes in the initial start codon; 
				if ($altStart == 1) {
					push(@cds_eff, [$cds_id, '3.7', "${refAA}${aa_pos}${altAA}"]); 
				} else {
					push(@cds_eff, [$cds_id, '3.1', "${refAA}${aa_pos}${altAA}"]); 
				}
				next; 
			}
			if ( $refAA eq '*' ) {
				# Changes in a reference stop codon; 
				if ( $altAA eq '*' ) {
					push(@cds_eff, [$cds_id, '3.6', "*${aa_pos}*"]); 
				} else {
					push(@cds_eff, [$cds_id, '3.2', "*${aa_pos}${altAA}"]); 
				}
				next; 
			} elsif ( $altAA eq '*' ) {
				# Change reference to stop codon; 
				push(@cds_eff, [$cds_id, '3.4', "${refAA}${aa_pos}${altAA}"]); 
				next; 
			}
			# Other changes in coding region; 
			if ( $refAA eq $altAA ) {
				push(@cds_eff, [$cds_id, '3.5', "${refAA}${aa_pos}${altAA}"]); 
			} else {
				push(@cds_eff, [$cds_id, '3.3', "${refAA}${aa_pos}${altAA}"]); 
			}
		}# End c2s position: for my $cp (@got_c2s)
		push(@back_eff, @cds_eff); 
	}
	@back_eff > 0 and return(@back_eff); # If the mutation relates to CDS, I don't care if it is in up/down-stream. 
	
	# For intron2scf : 
	my @intron_eff; 
	my @got_i2s = $mat_obj->switch_position(
		'qry2ref' => $agp_i2s_href, 
		'qryID'   => $chr_id, 
		'qryPos'  => $chr_pos, 
		'qryStr'  => '+'
	); 
	if (defined $got_i2s[0][0]) {
		# This means this position sits in at least one gene's upstream; 
		for my $cp (@got_i2s) {
			my ($t_id, $t_pos, $t_str) = @$cp; 
			$rb = $rb_ori; 
			$ab = $ab_ori; 
			my $t_id_ori=$t_id; 
			defined $in_intron{$t_id_ori}{'len'} or die "($t_id, $t_pos, $t_str) $chr_id $chr_pos |$in_intron{$t_id_ori}{'seq'}|\n"; 
			$t_id =~ s!_(\d+)$!!; 
			if ( $rb =~ m!^(InDel|Del|Ins)$!i or $ab =~ m!^(InDel|Del|Ins)$!i ) {
				push(@cds_eff, [$t_id, '2.2', "${t_id_ori}_${t_pos}_InDel"]); 
				next; 
			}
			if ($t_str eq '-') {
				# We need to get the reverse complement; 
				($rb) = &rev_comp($rb); 
				($ab) = &rev_comp($ab); 
			}
			if ( defined $opts{'chrFas'} ) {
				if ( $t_pos == 1 ) {
					my $refII = substr($in_intron{$t_id_ori}{'seq'}, 0, 2); 
					my $altII = $refII; 
					substr($refII, 0, 1) = $rb; 
					substr($altII, 0, 1) = $ab; 
					push(@intron_eff, [$t_id, '2.1', "${t_id_ori}_${refII}${t_pos}${altII}"]); 
				} elsif ( $t_pos == 2 ) {
					my $refII = substr($in_intron{$t_id_ori}{'seq'}, 0, 2); 
					my $altII = $refII; 
					substr($refII, 1, 1) = $rb; 
					substr($altII, 1, 1) = $ab; 
					push(@intron_eff, [$t_id, '2.1', "${t_id_ori}_${refII}${t_pos}${altII}"]); 
				} elsif ( $t_pos == $in_intron{$t_id_ori}{'len'}-1 ) {
					my $refII = substr($in_intron{$t_id_ori}{'seq'}, $in_intron{$t_id_ori}{'len'}-2, 2); 
					my $altII = $refII; 
					substr($refII, 0, 1) = $rb; 
					substr($altII, 0, 1) = $ab; 
					push(@intron_eff, [$t_id, '2.1', "${t_id_ori}_${refII}${t_pos}${altII}"]); 
				} elsif ( $t_pos == $in_intron{$t_id_ori}{'len'} ) {
					my $refII = substr($in_intron{$t_id_ori}{'seq'}, $in_intron{$t_id_ori}{'len'}-2, 2); 
					my $altII = $refII; 
					substr($refII, 1, 1) = $rb; 
					substr($altII, 1, 1) = $ab; 
					push(@intron_eff, [$t_id, '2.1', "${t_id_ori}_${refII}${t_pos}${altII}"]); 
				} else {
					push(@intron_eff, [$t_id, '2.2', "${t_id_ori}_$t_pos"]); 
				}
			} else {
				if ( $t_pos <= 2 or $t_pos >= $in_intron{$t_id_ori}{'len'}-1 ) {
					push(@intron_eff, [$t_id, '2.1', "${t_id_ori}_XX${t_pos}XX"]); 
				} else {
					push(@intron_eff, [$t_id, '2.2', "${t_id_ori}_$t_pos"]); 
				}
			}
		}
		push(@back_eff, @intron_eff); 
	}
	@back_eff > 0 and return(@back_eff); # If the mutation relates to intron, I don't care if it is in up/down-stream. 
	
	# For up2scf : Check if position is in upstream; 
	my @up_eff; 
	my @got_u2s = $mat_obj->switch_position(
		'qry2ref' => $agp_u2s_href, 
		'qryID'   => $chr_id, 
		'qryPos'  => $chr_pos, 
		'qryStr'  => '+'
	); 
	if (defined $got_u2s[0][0]) {
		# This means this position sits in at least one gene's upstream; 
		for my $cp (@got_u2s) {
			my ($t_id, $t_pos, $t_str) = @$cp; 
			push(@up_eff, [$t_id, '1.1', $t_pos]); 
		}
		push(@back_eff, @up_eff); 
	}
	# For down2scf : Check if position is in downstream; 
	my @down_eff; 
	my @got_d2s = $mat_obj->switch_position(
		'qry2ref' => $agp_d2s_href, 
		'qryID'   => $chr_id, 
		'qryPos'  => $chr_pos, 
		'qryStr'  => '+'
	); 
	if (defined $got_d2s[0][0]) {
		# This means this position sits in at least one gene's upstream; 
		for my $cp (@got_d2s) {
			my ($t_id, $t_pos, $t_str) = @$cp; 
			push(@down_eff, [$t_id, '1.2', $t_pos]); 
		}
		push(@back_eff, @down_eff); 
	}
	# If there is no effect found, this SNP locates in intergenic region; 
	if (@back_eff == 0) {
		push(@back_eff, ['NULL', '1.3']); 
	}
	return(@back_eff); 
}# snp_eff() 

sub rev_comp {
	my @back; 
	for (@_) {
	defined $_ or die "aaa\n"; 
		my $ts = reverse($_); 
		$ts =~ tr/acgturykmbvdhACGTURYKMBVDHwWsSnN/tgcaayrmkvbhdTGCAAYRMKVBHDwWsSnN/; 
		push(@back, $ts); 
	}
	return(@back); 
}# rev_comp() 
##### Get 

sub load_gff {
	my ($iFh, $gff_href, $seq_href) = @_; 
	$gff_href //= {}; 
	$seq_href //= {}; 
	my %in_gff = %$gff_href; 
	my %in_seq = %$seq_href; 
	my ( $in_gff_href, $in_seq_href ) = $gff_obj->read_gff3File('gffFH'=>$iFh, 'saveFa'=>0, 'top_hier'=>$glob{'gff_top_hier'}); 
	close($iFh); # 
	%in_gff = %$in_gff_href; 
	%in_seq = %$in_seq_href; 
	%$gff_href = %in_gff; 
	%$seq_href = %in_seq; 
	return($gff_href, $seq_href); 
}# load_gff() 

sub retrive_cds {
	my %str2num = qw(
	 +      1
	 -     -1
	 1      1
	 -1    -1
	 0     -1
	 plus   1
	 minus -1
	); 
	TOPID:
	for my $topID ( sort { $in_gff{'ID2lineN'}{$a} <=> $in_gff{'ID2lineN'}{$b} } keys %{$in_gff{'lineN_group'}} ) {
		my @posi_top; 
		my @posi_cds; 
		my $curLnNum = $in_gff{'lineN_group'}{$topID}{'curLn'}[0]; 
		my $top_str = ''; 
		my $top_chr = ''; 
		my $top_name = ''; 
		my %curLnH; 
		if ( defined $in_gff{'lineN2hash'}{$curLnNum} ) {
			%curLnH = %{ $in_gff{'lineN2hash'}{$curLnNum} }; 
			$top_str = $str2num{ lc($curLnH{'strand'}) } // ''; 
			$top_chr = $curLnH{'seqID'} // ''; 
			$top_name = $curLnH{'attrib'}{'featID'} // ''; 
			@posi_top = ( [$curLnH{'start'}, $curLnH{'end'}] ); 
		}
		OFFID:
		for my $offLnNum ( @{$in_gff{'lineN_group'}{$topID}{'offLn'}} ) {
			defined $in_gff{'lineN2hash'}{$offLnNum} or next OFFID; 
			my %offLnH = %{ $in_gff{'lineN2hash'}{$offLnNum} }; 
			$offLnH{'type'} =~ m!^cds$!i or next OFFID; 
			$top_str eq '' and $top_str =  $str2num{ lc($offLnH{'strand'}) } // ''; 
			 $top_chr eq '' and $top_chr = $offLnH{'seqID'} // ''; 
			 if ( keys %{$offLnH{'attrib'}{'parentID'}} > 0 ) {
				my @ta1 = sort keys %{ $offLnH{'attrib'}{'parentID'} };
				my $ta1_txt = join(";", @ta1);
				if ( $top_name eq '' ) {
					$top_name = $ta1_txt; 
				} else {
					$top_name eq $ta1_txt or &stopErr("[Err] Unequal top_name : [ $top_name, $ta1_txt ]\n"); 
				}
			 }
			 my @ta = split(/\t/, $in_gff{'lineN2line'}{$offLnNum}); 
			 my $t_frame = $ta[7]; 
			 $t_frame eq '.' and $t_frame = 0; 
			 $t_frame ++; 
			 push( @posi_cds, [ $offLnH{'start'}, $offLnH{'end'}, $t_frame] ); 
		}# End for my $offLnNum
		@posi_cds > 0 or next TOPID; 
		$top_str eq '' and do { &tsmsg("[Wrn] No strand information for topID=[$topID]\n"); $top_str = 1; };
		$top_chr eq '' and do { &stopErr("[Err] No top_chr found for [$topID]\n"); };
		defined $in_seq{$top_chr} or do { &tsmsg("[Wrn] No [$top_chr] sequence found for [$topID]\n"); next TOPID; };
		# Check @posi_cds
		for my $tr (@posi_cds) {
			$tr->[0] > $tr->[1] and &stopErr("[Err] I can't accept start[$tr->[0]] > end[$tr->[1]] in gff3 file.\n"); 
		}
		if ($top_str == -1) {
			# Be careful when the CDS loci should not ordered by position; 
			@posi_cds = sort { $b->[0] <=> $a->[0] } @posi_cds;
		}
		
		# Generate AGP file for CDS; 
		my %h1; # This is used for cds2scf.agp; 
		$h1{'eleN'} = 0; 
		for my $tr (@posi_cds) {
			$h1{'eleN'} ++; 
			$h1{'prevE'} //= 0; 
			print {$glob{'fh_cds2Scf_agp'}} join("\t", 
				$top_name, 
				$h1{'prevE'}+1, 
				$h1{'prevE'}+$tr->[1]-$tr->[0]+1,
				$h1{'eleN'},
				'W', 
				$top_chr, 
				$tr->[0], 
				$tr->[1], 
				( $top_str == -1 ) ? '-' : '+'
			)."\n"; 
			if ( $h1{'eleN'} > 1 ) {
				# For intron2Scf
				my $iname = $top_name . "_" . ($h1{'eleN'}-1); 
				if ( $top_str == -1 ) {
					print {$glob{'fh_intron2Scf_agp'}} join("\t", 
						$iname, 
						1, 
						$h1{'intronS'}-$tr->[1], 
						1, 
						'W', 
						$top_chr, 
						$tr->[1]+1, 
						$h1{'intronS'}, 
						'-'
					)."\n"; 
					$in_intron{$iname}{'len'} = $h1{'intronS'}-$tr->[1]; 
					if (defined $opts{'chrFas'}) {
						my $iseq = substr($in_seq{$top_chr}{'seq'}, $tr->[1], $h1{'intronS'}-$tr->[1]); 
						($iseq) = &rev_comp($iseq); 
						$in_intron{$iname}{'seq'} = $iseq; 
					}
				} else {
					print {$glob{'fh_intron2Scf_agp'}} join("\t", 
						$top_name . "_" . ($h1{'eleN'}-1), 
						1, 
						$tr->[0]-$h1{'intronS'}, 
						1, 
						'W', 
						$top_chr, 
						$h1{'intronS'}, 
						$tr->[0]-1, 
						'+'
					)."\n"; 
					$in_intron{$iname}{'len'} = $tr->[0]-$h1{'intronS'}; 
					if (defined $opts{'chrFas'}) {
						my $iseq = substr($in_seq{$top_chr}{'seq'}, $h1{'intronS'}-1, $tr->[0]-$h1{'intronS'}); 
						$in_intron{$iname}{'seq'} = $iseq; 
					}
				}
			}
			$h1{'prevE'} = $h1{'prevE'}+$tr->[1]-$tr->[0]+1; 
			$h1{'intronS'} = ( $top_str == -1 ) ? $tr->[0]-1 : $tr->[1]+1 ; 
		}# End for my $tr (@posi_cds) 
		if ($top_str == -1) {
			print {$glob{'fh_up2Scf_agp'}} join("\t", 
				$top_name, 
				1, 
				$glob{'upstreamLen'}, 
				1, 
				'W', 
				$top_chr, 
				$posi_cds[0][1]+1, 
				$posi_cds[0][1]+$glob{'upstreamLen'}, 
				'+'
			)."\n"; 
			print {$glob{'fh_down2Scf_agp'}} join("\t", 
				$top_name, 
				1, 
				$glob{'downstreamLen'}, 
				1, 
				'W', 
				$top_chr, 
				$posi_cds[-1][0]-$glob{'downstreamLen'}, 
				$posi_cds[-1][0]-1, 
				'-'
			)."\n"; 
		} else {
			print {$glob{'fh_up2Scf_agp'}} join("\t", 
				$top_name, 
				1, 
				$glob{'upstreamLen'}, 
				1, 
				'W', 
				$top_chr, 
				$posi_cds[0][0]-$glob{'upstreamLen'}, 
				$posi_cds[0][0]-1, 
				'-'
			)."\n"; 
			print {$glob{'fh_down2Scf_agp'}} join("\t", 
				$top_name, 
				1, 
				$glob{'downstreamLen'}, 
				1, 
				'W', 
				$top_chr, 
				$posi_cds[-1][1]+1, 
				$posi_cds[-1][1]+$glob{'downstreamLen'}, 
				'+'
			)."\n"; 
		}
		
		# Generate CDS sequences and frames for further usage; 
		if ( defined $opts{'chrFas'} ) {
			# This means that %in_seq stores the scf sequences; 
			my @sub_seqs; 
			for my $tr (@posi_cds) {
				push( @sub_seqs, substr($in_seq{$top_chr}{'seq'}, $tr->[0]-1, $tr->[1]-$tr->[0]+1) );
			}
			$top_str == -1 and @sub_seqs = &rev_comp(@sub_seqs);
			my $final_seq = join('', @sub_seqs); 
			$in_cds{$top_name}{'seq'} = $final_seq; 
			$in_cds{$top_name}{'frame'} = $posi_cds[0][2]; 
		}
	}# End for my $topID 
	close($glob{'fh_cds2Scf_agp'}); 
	close($glob{'fh_intron2Scf_agp'}); 
	close($glob{'fh_up2Scf_agp'}); 
	close($glob{'fh_down2Scf_agp'}); 
	%agp_cds2Scf  = %{ &fileSunhh::load_agpFile( $glob{'file_cds2Scf_agp'} ) }; 
	%agp_intron2Scf  = %{ &fileSunhh::load_agpFile( $glob{'file_intron2Scf_agp'} ) }; 
	%agp_up2Scf   = %{ &fileSunhh::load_agpFile( $glob{'file_up2Scf_agp'} ) }; 
	%agp_down2Scf = %{ &fileSunhh::load_agpFile( $glob{'file_down2Scf_agp'} ) }; 
	return(); 
}# retrive_cds() 


