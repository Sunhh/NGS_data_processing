#!/usr/bin/perl -w
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 

use Parallel::ForkManager;


use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
 "in_gff:s", # in_gff is similar to mcscan's .bed file, with format as : (OO)N+ \t GeneID \t Start \t End 
 "in_aln:s", # This is the output ".collinearity" file of mcscanX, which should be similar to .align of mcscan or previous mcscanX version. 
 "in_blast:s", # This is the input ".blast" file of mcscanX, 
 "out:s", # output filename 
 "aln2list!", "addChr!", "tgt_gff:s", "srt_by:s", 
 "aln2table!", 
 "glist2pairs!", "in_glist:s", "pivot_pat:s", "target_pat:s", 
 "slctBlk!", "slct_list:s", "slct_colN:i", "slct_type:s", 
 "glist2html!", "pivot_chrID:s", 
 "classDupGene!", "max_eval:f", "min_cscore:f", "max_proxN:i", "max_tandN:i", "tand_list:s", 
 "filterBlast!", "rm_repeat!", "repeat_eval:f", "rm_tandem!", 
 "add_KaKs!", "in_pair_list:s", "fas_cds:s", "fas_prot:s", "out_pref:s","alnMethod:s", "ncpu:i", 
 "help!", 
); 

sub usage {
	print <<HH;
####################################################################################################
# perl $0 
# 
# -help 
#####  Basic input files 
# -in_gff     similar to mcscan's .bed file, with format as : (OO)N+ \\t GeneID \\t Start \\t End
#               Here 'OO' is a two-characters organism name, 'N+' is the chromosome number. 
# -in_aln     the output ".collinearity" file of mcscanX, which should be similar to .align of mcscan or previous mcscanX version. 
# -in_blast   the input ".blast" file of mcscanX, formatted as "blast -m8", 
#               "gen1 \\t gen2 \\t score \\t blkLen \\t mismatch \\t gapopen \\t qS \\t qE \\t sS \\t sE \\t evalue \\t bitscore"
#####  
# -out        [outfilename] Default output to STDOUT. 
#####  Functions 
# -aln2list       [Boolean] Extract a tabular list from -in_aln according to -in_aln. 
#                   Need -in_gff , -in_aln 
#                   This list can be viewed in excel, and each column contain an aligned block. 
#   -addChr       [Boolean] Tell if add chromID to block genes. 
#   -tgt_gff      [target.gff] Provide source of blocks. 
#   -srt_by       ['min'] Sort the gene blocks by 'min|Q1|Q3|interval_mean|COUNT'
# 
# -aln2table      [Boolean] Reformat aligned blocks into one line. 
#                   Need -in_gff , -in_aln 
#                   Headers: BlkID / Chrom1 / Start1 / End1 / Chrom2 / Start2 / End2 / Strand / AlnScore / AlnEvalue / AlnNumber / Gene1 / Gene2 / 
#                            Ka / Ks / KaKs / AVG_Ks / Med_Ks / UsedKs / AVG_Ka / Med_Ka / UsedKa / AVG_KaKs / Med_KaKs / UsedKaKs
# 
# -glist2pairs    [Boolean] Extract non-redundant gene pairs in gene list. 
#                   Here the gene list is the output of -aln2list. 
#                   Please note that this pairs' set might be larger than pairs' set directly shown in .aln/.collinearity file, 
#                   because if there are only A-B and A-C relationship in .collinearity file, we will get B-C pair in addition. 
#                   Need -in_glist 
#   -in_glist     [filename] Must be given if '-glist2pairs' assigned. 
#   -pivot_pat    [pattern_string] ''; Only check lines with 'pivot_pat' pattern in the 1st-column ('Chromosome', colN=0); 
#   -target_pat   [pattern_string] ''; Only check target genes with 'target_pat' pattern in the target gene columns ('Blocks', colN >= 6)
# 
# -slctBlk        [Boolean] Extract alignment-block according to '-slct_list' and '-slct_type'
#                   Need -in_aln . 
#   -slct_list    [filename] Must be given if '-slctBlk' assigned. 
#   -slct_colN    [0] 0-based selected column number. 
#   -slct_type    ['blkID'] Type of selected column. Could be 'blkID|chrID'
#
# -glist2html     [Boolean] Convert -aln2list output to .html format
#                   Need -in_glist . 
#   -in_glist     [filename] Must be given if '-glist2html' assigned. 
#   -pivot_chrID  [ChromID] The chromosome ID which is used as .html file reference. Must be given if '-glist2html' assigned. 
# 
# -classDupGene   [Boolean] Classify duplicated genes into classes 'Singleton|Dispersed|Proximal|Tandem|WGD or segmental'
#                   Singleton        : [0] 
#                   Dispersed        : [1] 
#                   Proximal         : [2] 
#                   Tandem           : [3] 
#                   WGD or segmental : [4] 
#                   Need -in_gff , -in_aln , -in_blast . 
#   -max_eval     [-1] Filter blast gene pairs by e-value. Only pairs with e-value <= \$max_eval is accepted. -1 means no filter. 
#   -min_cscore   [-1] Filter blast gene pairs by cscore. cscore(A,B) = score(A,B) / max(best score for A, best score for B)
#                   http://www.sciencemag.org/cgi/content/abstract/317/5834/86
#                   https://github.com/tanghaibao/quota-alignment/blob/master/scripts/blast_to_raw.py
#   -max_proxN    [20] The max gene number allowed between proximal duplicated genes. 
#   -max_tandN    [2] The max gene number allowed between tandem duplicated genes. 
#   -tand_list    [outFileName] If given, write a table to 'outFileName' with foramt: 
#                   Grp0_gen1 \\t Grp0_gen2 \\t Grp0_gen3 ...
#                   Grp1_gen1 \\t Grp1_gen2 \\t Grp1_gen3 ... 
#
# -filterBlast    [Boolean] Filter .blast file. 
#                   Need -in_blast 
#   -rm_repeat    [Boolean]
#    -repeat_eval [0.05] 
#   -rm_tandem    [Boolean] Should be given with -in_gff 
#   -max_tandN
#   -min_cscore
#   -max_eval
#
# -add_KaKs       [Boolean] Invoke tanghaibao-jcvi python to do this. 
#                   Need -in_aln / -in_pair_list 
#   -in_pair_list [filename] Format: 'gene1 \\t gene2'
#   -fas_cds      [filename] fasta format of cds file. 
#   -fas_prot     [filename] fasta format of protein file. 
#   -out_pref     [Prefix] Output prefix. Default is 'paired'. 
#   -alnMethod    ['muscle'] Could also be 'clustalw'
#   -ncpu         [1] Use multiple-CPUs to compute KaKs, because it is too slow. 
####################################################################################################
HH
	exit 1; 
}

$opts{'help'} and &usage(); 
(keys %opts) == 0 and &usage(); 

####################################################################################################
# Set basic parameters. 
####################################################################################################
my $outFh = \*STDOUT; 
defined $opts{'out'} and $outFh = &openFH($opts{'out'}, '>'); 
$opts{'srt_by'} //= 'min'; 


####################################################################################################
# Call different sub-funcitons. 
####################################################################################################

if ( $opts{'aln2list'} ) {
	&msc_aln2list(
	 'in_gff'=>$opts{'in_gff'}, 
	 'in_aln'=>$opts{'in_aln'}, 
	 'addChr'=>$opts{'addChr'}, 
	 'tgt_gff'=>$opts{'tgt_gff'}, 
	 'srt_by'=>$opts{'srt_by'}
	); 
} elsif ( $opts{'aln2table'} ) {
	&msc_aln2table(
	 'in_gff'=>$opts{'in_gff'}, 
	 'in_aln'=>$opts{'in_aln'}, 
	); 
} elsif ( $opts{'glist2pairs'} ) {
	&msc_glist2pairs(
	 'in_glist'=>$opts{'in_glist'}, 
	 'pivot_pat'=>$opts{'pivot_pat'}, 
	 'target_pat'=>$opts{'target_pat'}, 
	); 
} elsif ( $opts{'slctBlk'} ) {
	&mcs_blkByList(
	 'in_aln'=>$opts{'in_aln'}, 
	 'slct_list'=>$opts{'slct_list'}, 
	 'slct_colN'=>$opts{'slct_colN'}, 
	 'slct_type'=>$opts{'slct_type'}, 
	); 
} elsif ( $opts{'glist2html'} ) {
	&mcs_glist2html(
	 'in_glist'=>$opts{'in_glist'}, 
	 'pivot_chrID'=>$opts{'pivot_chrID'}, 
	); 
} elsif ( $opts{'classDupGene'} ) { 
	&mcs_classDup(
	 'in_gff'=>$opts{'in_gff'}, 
	 'in_aln'=>$opts{'in_aln'}, 
	 'in_blast'=>$opts{'in_blast'}, 
	 'max_eval'=>$opts{'max_eval'}, 
	 'min_cscore'=>$opts{'min_cscore'}, 
	 'max_proxN'=>$opts{'max_proxN'}, 
	 'max_tandN'=>$opts{'max_tandN'}, 
	 'tand_list'=>$opts{'tand_list'}, 
	); 
} elsif ( $opts{'filterBlast'} ) { 
	&mcs_filterBlast(
	 'in_gff'=>$opts{'in_gff'}, 
	 'in_blast'=>$opts{'in_blast'}, 
	 'max_eval'=>$opts{'max_eval'}, 
	 'min_cscore'=>$opts{'min_cscore'}, 
	 'rm_repeat'=>$opts{'rm_repeat'}, 
	 'repeat_eval'=>$opts{'repeat_eval'}, 
	 'rm_tandem'=>$opts{'rm_tandem'}, 
	 'max_tandN'=>$opts{'max_tandN'}, 
	)
} elsif ( $opts{'add_KaKs'} ) {
	$opts{'out_pref'} //= 'paired'; 
	&mcs_addKaKs(
	 'in_aln'=>$opts{'in_aln'}, 
	 'in_pair_list'=>$opts{'in_pair_list'}, 
	 'fas_cds'=>$opts{'fas_cds'}, 
	 'fas_prot'=>$opts{'fas_prot'}, 
	 'out_pref'=>$opts{'out_pref'}, 
	 'alnMethod'=>$opts{'alnMethod'}, 
	 'ncpu' => $opts{'ncpu'}, 
	); 
} else {
	&usage(); 
}
####################################################################################################
# Sub-functions to be called 
####################################################################################################

sub mcs_addKaKs {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	$parm{'fas_prot'} //= ''; 
	$parm{'out_pref'} //= 'paired'; 
	$parm{'ncpu'} //= 1; 
	defined $parm{'fas_cds'} or &stopErr("[Err] -fas_cds must be given.\n"); 
	
	if ( defined $parm{'in_pair_list'} ) {
		if ( $parm{'ncpu'} > 1 ) {
			my @sub_dirs; 
			my $ipFh = &openFH($parm{'in_pair_list'}, '<'); 
			my @in_lines = <$ipFh>; 
			close($ipFh); 
			my $lnPerFile = scalar(@in_lines) / $parm{'ncpu'} ; 
			$lnPerFile = ( $lnPerFile > int($lnPerFile) ) ? $lnPerFile+1 : $lnPerFile ; 
			for (my $i=0; $i<@in_lines; $i+=$lnPerFile) {
				my $e = $i+$lnPerFile-1; $e>$#in_lines and $e = $#in_lines; 
				my $tmpDir = &fileSunhh::new_tmp_dir(); 
				push(@sub_dirs, $tmpDir); 
				mkdir($tmpDir); 
				my $opFh = &openFH("$tmpDir/pair_lis", '>'); 
				for my $j ($i .. $e) {
					print {$opFh} $in_lines[$j]; 
				}
				close($opFh); 
			}
			my $pm = new Parallel::ForkManager($parm{'ncpu'}); 
			for my $sd ( @sub_dirs ) {
				my $pid = $pm->start and next; 
				my %sub_parm = %parm; 
				for my $tk (qw/in_aln in_pair_list fas_cds fas_prot/) {
					defined $sub_parm{$tk} or next; 
					$sub_parm{$tk} = fileSunhh::_abs_path( $sub_parm{$tk} ); 
				}
				chdir($sd); 
				$sub_parm{'in_pair_list'} = 'pair_lis'; 
				$sub_parm{'out_pref'} = 'paired'; 
				&_lis2ks(%sub_parm); 
				$pm->finish; 
			}
			$pm->wait_all_children; 
			my $oksFh = &openFH("$parm{'out_pref'}.ks", '>'); 
			my $has_head = 0; 
			for my $sd ( @sub_dirs ) {
				my $sksFh = &openFH("$sd/paired.ks", '<'); 
				while (<$sksFh>) { $. == 1 and $has_head == 1 and next; $has_head = 1; chomp; print {$oksFh} "$_\n"; } 
				close($sksFh); 
			}
			close($oksFh); 
			&fileSunhh::_rmtree(\@sub_dirs); 
		} else {
			&_lis2ks(%parm); 
		}
	} elsif ( $parm{'in_aln'} ) {
		my ($alnInfo) = &_readInAln( $parm{'in_aln'} ); 
		$parm{'in_pair_list'} = "$parm{'out_pref'}_pair.lis"; 
		my $lisFh = &openFH($parm{'in_pair_list'}, '>'); 
		for my $ar1 (@$alnInfo) {
			for my $ar2 (@{$ar1->{'pair'}}) {
				my ($g1, $g2) = ($ar2->[0], $ar2->[1]); 
				print {$lisFh} "$g1\t$g2\n"; 
			}
		}
		close($lisFh); 
		my ($header, $pair2ks); 
		if ( $parm{'ncpu'} > 1 ) {
			my %tmp_parm = %parm; 
			&mcs_addKaKs(%tmp_parm); 
			( $header, $pair2ks ) = &_readInKsTab( "$tmp_parm{'out_pref'}.ks" ); 
		} else {
			( $header, $pair2ks ) = &_lis2ks(%parm); 
		}
		my $oksFh = &openFH("$parm{'in_aln'}.ks", '>'); 
		for my $ar1 (@$alnInfo) {
			for my $ln (split(/\n/, $ar1->{'text'})) {
				chomp($ln); 
				if ( $ln =~ m/^\s*(#|$)/ ) {
					print {$oksFh} "$ln\n"; 
					next; 
				}
				my @ta = split(/\t/, $ln); 
				print {$oksFh} join( "\t", @ta[0..3], @{ $pair2ks->{$ta[1]}{$ta[2]} } )."\n"; 
			}
		}
		close ($oksFh); 
	}
	return ; 
}# mcs_addKaKs() 

# Input file format : 
#   gene1   gene2   yn_ks   yn_ka   ng_ks   ng_ka
#   Cma_000007      Cma_000973      0.5258  0.0154  0.5387  0.0148
#   Cma_000009      Cma_000975      0.2844  0.0044  0.2966  0.0043
sub _readInKsTab {
	my $inFh = &openFH(shift, '<'); 
	my (@header, %pair2ks); 
	while (<$inFh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		if ($ta[0] eq 'gene1') {
			push(@header, @ta); 
			next; 
		}
		($ta[0] ne '' and $ta[1] ne '') or do { &tsmsg("[Err] bad line in .ks.tab file: $_\n"); next; }; 
		$pair2ks{$ta[0]}{$ta[1]} = [ @ta[2..$#ta] ]; 
	}
	close($inFh); 
	return(\@header, \%pair2ks); 
}# _readInKsTab () 


# Input file format : 
#   name,yn_ks,yn_ka,ng_ks,ng_ka
#   Cla007160;Cla010603,0.4049,0.2933,0.4497,0.2698
sub _readInKsLis {
	my $inFh = &openFH(shift, '<'); 
	my (@header, %pair2ks); 
	while (<$inFh>) {
		chomp; 
		my @ta = split(/,/, $_); 
		if ($ta[0] eq 'name') {
			push(@header, qw/gene1 gene2/, @ta[1..$#ta]); 
			next; 
		}
		$ta[0] =~ m!^([^\s;]+);([^\s;]+)$! or &stopErr("[Err] bad line in .ks file: $_\n"); 
		$pair2ks{$1}{$2} = [ @ta[1..$#ta] ]; 
	}
	close($inFh); 
	return(\@header, \%pair2ks); 
}# _readInKsLis () 

sub _lis2ks {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	$parm{'fas_prot'} //= ''; 
	$parm{'out_pref'} //= 'paired'; 
	$parm{'alnMethod'} //= 'muscle'; 
	defined $parm{'fas_cds'} or &stopErr("[Err] -fas_cds must be given.\n"); 
	defined $parm{'in_pair_list'} or &stopErr("[Err] -in_pair_list is needed.\n"); 
	
	&exeCmd_1cmd("python -m jcvi.apps.ks prepare $parm{'in_pair_list'} $parm{'fas_cds'} $parm{'fas_prot'} -o $parm{'out_pref'}.cds.fasta") and &stopErr("[Err] Failed prepare\n"); 
	my $med_prot = ( $parm{'fas_prot'} eq '' ) ? "" : "$parm{'out_pref'}.cds.fasta.pep" ; 
	&exeCmd_1cmd("python -m jcvi.apps.ks calc    --msa=$parm{'alnMethod'} $med_prot $parm{'out_pref'}.cds.fasta -o $parm{'out_pref'}.cds.fasta.ks"); 
	my ( $header, $pair2ks ) = &_readInKsLis( "$parm{'out_pref'}.cds.fasta.ks" ); 
	my $paFh = &openFH($parm{'in_pair_list'}, '<'); 
	my $oksFh = &openFH("$parm{'out_pref'}.ks", '>'); 
	# print {$oksFh} join("\t", qw/gene1 gene2/, @$header[1 .. $#{$header}])."\n"; 
	print {$oksFh} join("\t", @$header[0 .. $#{$header}])."\n"; 
	while (<$paFh>) {
		m/^\s*(#|$)/ and next; 
		chomp; s!\s+$!!; 
		my @ta = split(/\t/, $_); 
		my ($g1, $g2) = @ta[0,1]; 
		defined $pair2ks->{$g1}{$g2} or do { &tsmsg("[Err] No Ks result found for line: $_\n"); next; }; 
		print {$oksFh} join("\t", $g1, $g2, @{$pair2ks->{$g1}{$g2}})."\n"; 
	}
	close ($paFh); 
	
	return ($header, $pair2ks);
}# _lis2ks() 


sub mcs_filterBlast {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	$parm{'max_eval'}    //= -1; 
	$parm{'min_cscore'}  //= -1; 
	$parm{'rm_repeat'}   //= undef(); 
	$parm{'rm_tandem'}   //= undef(); 
	$parm{'repeat_eval'} //= 0.05; 
	$parm{'max_tandN'}   //= 2; 
	my ($bn6Info) = &_readInBn6( $parm{'in_blast'} ); 
	my ($chr2gen, $gen2loc); 
	($chr2gen, $gen2loc) = &_readInGff($parm{'in_gff'}) if ( defined $parm{'in_gff'} and $parm{'in_gff'} ne '' ) ; 
	
	$bn6Info = &_filter_eval( $bn6Info, $parm{'max_eval'} ) if ( $parm{'max_eval'} >= 0 ); 
	$bn6Info = &_filter_tand( $bn6Info, $parm{'max_tandN'}, $gen2loc ) if ( $parm{'rm_tandem'} ); 
	$bn6Info = &_filter_repeat( $bn6Info, $parm{'repeat_eval'} ) if ( $parm{'rm_repeat'} ); 
	$bn6Info = &_filter_cscore( $bn6Info, $parm{'min_cscore'}, $gen2loc ) if ( $parm{'min_cscore'} > 0 ); 
	
	for my $ar1 (@$bn6Info) {
		print {$outFh} join("\t", @{$ar1->{'arr'}})."\n"; 
	}
	
	return; 
}# mcs_filterBlast() 

sub mcs_classDup {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	$parm{'max_eval'}   //= -1; 
	$parm{'min_cscore'} //= -1; 
	$parm{'max_proxN'}  //= 20; 
	$parm{'max_tandN'}  //= 2; 
	
	my ($chr2gen, $gen2loc) = &_readInGff($parm{'in_gff'}); 
	my %gidx = %{ &_index_glist($chr2gen) }; 
	
	my ($alnInfo) = &_readInAln( $parm{'in_aln'} ); 
	my ($bn6Info) = &_readInBn6( $parm{'in_blast'} ); 
	$bn6Info = &_filter_eval( $bn6Info, $parm{'max_eval'} ) if ($parm{'max_eval'} >= 0); 
	$bn6Info = &_filter_cscore( $bn6Info, $parm{'min_cscore'}, $gen2loc ) if ($parm{'min_cscore'} > 0) ; 
	my ($id_grp) = &_group_id( $bn6Info ); 
	# [[grp0_gen1, grp0_gen2, ...], [grp1_gen1, grp2_gen2, ...], ...]
	my $prox_gene = &_tandemGrp(
	 'id_grp'=>$id_grp, 
	 'max_distN'=>$parm{'max_proxN'}, 
	 'gen2loc'=>$gen2loc, 
	 'gindex'=>\%gidx, 
	); 
	my $tand_gene = &_tandemGrp(
	 'id_grp'=>$id_grp, 
	 'max_distN'=>$parm{'max_tandN'}, 
	 'gen2loc'=>$gen2loc, 
	 'gindex'=>\%gidx, 
	); 
	
	if ( defined $parm{'tand_list'} ) {
		my $otFh = &openFH($parm{'tand_list'}, '>'); 
		for my $ar1 (sort { $gidx{$a->[0]} <=> $gidx{$b->[0]} } @$tand_gene) {
			print {$otFh} join("\t", @$ar1)."\n"; 
		}
		close($otFh); 
	}
	
	my %gene_tag; 
	# {gene_ID}{0/1/2/3/4} => 1,
	# 0-singletons,
	# 1-dispersed duplicates,
	# 2-proximal duplicates,
	# 3-tandem duplicates,
	# 4-segmental/WGD duplicates.
	
	for my $genID ( keys %$gen2loc ) { $gene_tag{$genID}{0} = 0; } 
	for my $ar1 ( @$bn6Info ) {
		$gene_tag{ $ar1->{'arr'}[0] }{1} = 1; 
		$gene_tag{ $ar1->{'arr'}[1] }{1} = 1; 
	}
	for my $ar1 ( @$prox_gene ) {
		for my $genID (@$ar1) { $gene_tag{$genID}{2} = 1; }
	}
	for my $ar1 ( @$tand_gene ) {
		for my $genID (@$ar1) { $gene_tag{$genID}{3} = 1; }
	}
	for my $ar1 ( @$alnInfo ) {
		for my $ar2 ( @{$ar1->{'pair'}} ) {
			$gene_tag{ $ar2->[0] }{4} = 1; 
			$gene_tag{ $ar2->[1] }{4} = 1; 
		}
	}
	
	print {$outFh} join("\t", qw/GeneID ClassN C0_1_2_3_4/)."\n"; 
	for my $genID (sort { $gidx{$a}<=>$gidx{$b} } keys %$gen2loc) {
		my %tt_hash = %{ $gene_tag{ $genID } }; 
		my $final_tag = 0; 
		map { $final_tag < $_ and $final_tag = $_; } grep { defined $tt_hash{ $_ } and $tt_hash{$_} > 0 } keys %tt_hash ; 
		for (qw/0 1 2 3 4/) { $tt_hash{$_} //= 0; }
		print {$outFh} join( "\t", $genID, $final_tag, join("_", @tt_hash{qw/0 1 2 3 4/}) )."\n"; 
	}
}# mcs_classDup() 

sub mcs_glist2html {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	defined $parm{'pivot_chrID'} or &stopErr("[Err] -pivot_chrID not defined.\n"); 
	my $chrID = $parm{'pivot_chrID'}; 
	
	my ($glist_arr) = &_readInGlist($parm{'in_glist'}); # 
	
	my $mid_head = <<MIDHEAD; 
<tr align='center'>
MIDHEAD
	my $nbsp = '&nbsp;'; 
	my %trans_code; 
	$trans_code{'||'} = &_write_td("|$nbsp|"); 
	$trans_code{'null'} = &_write_td("$nbsp"); 
	$trans_code{'space'} = &_write_td("$nbsp$nbsp"); 
	
	my %color; 
	$color{'pivot'}  = '#dddddd';
	$color{'tandem'} = '#ff0000';
	$color{'paired'} = '#ffff99';
	$color{'start'}  = '#ff00ff';
	$color{'end'}    = '#0000ff';

	print {$outFh} <<HEAD; 
<html>
<table cellspacing='0' cellpadding='0' align='left'>
<tr align='center'>
<td>Intra dup depth</td>
$trans_code{space}
<td>Inter dup depth</td>
$trans_code{space}
<td>Inter dup taxes</td>
$trans_code{space}
<td>Start Position</td>
$trans_code{space}
<td>Pivot gene:$chrID</td>
$trans_code{space}
<td align='left' colspan='152'>Collinear blocks</td>
</tr>
HEAD
	
	for my $ar1 ( @{$glist_arr} ) {
		$ar1->[0] eq $chrID or next; 
		print {$outFh} "<tr align='center'>"; 
		print {$outFh} "<td>$ar1->[1]</td>"; # intra-depth
		print {$outFh} "$trans_code{space}<td align='right'>$ar1->[2]</td>"; # inter-depth
		$ar1->[3] eq '' and $ar1->[3] = ${nbsp}; # No inter-tax exists, use space ' ' to replace NULL. 
		print {$outFh} "$trans_code{space}<td align='right'>$ar1->[3]</td>"; # inter-taxes
		print {$outFh} "$trans_code{space}<td align='right'>$ar1->[5]</td>"; # Start-Position
		print {$outFh} "$trans_code{space}<td bgcolor='$color{pivot}'>$ar1->[4]</td>"; # Pivot-GeneID 
		for (my $i=6; $i<@$ar1; $i++) {
			if ( $ar1->[$i] eq '' ) {
				print {$outFh} "$trans_code{space}<td>${nbsp}</td>"; 
			} elsif ( $ar1->[$i] eq '||' ) {
				print {$outFh} "$trans_code{space}<td>|${nbsp}|</td>"; 
			} else {
				my $tag = 'paired'; 
				if ( $ar1->[$i] =~ s!\-e$!! ) {
					$tag = 'end'; 
				} elsif ( $ar1->[$i] =~ s!\-s$!! ) {
					$tag = 'start'; 
				}
				print {$outFh} "$trans_code{space}<td bgcolor='$color{$tag}'>$ar1->[$i]</td>"; 
			}
		}
		print {$outFh} "$trans_code{space}<td>${nbsp}</td>"; 
		print {$outFh} "</tr>\n"; 
	}
}# mcs_glist2html() 

sub mcs_blkByList {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	
	my ($alnInfo) = &_readInAln($parm{'in_aln'}); 
	my $inFh = &openFH($parm{'slct_list'}, '<'); 
	$parm{'slct_colN'} //= 0; 
	$parm{'slct_type'} //= 'blkID'; 
	my %is_want; 
	while (<$inFh>) {
		chomp; m/^\s*$/ and next; 
		my @ta = split(/\t/, $_); 
		$is_want{ $ta[ $parm{'slct_colN'} ] } = 1; 
	}
	close($inFh); 
	defined $alnInfo->[0]{'text'} and print {$outFh} $alnInfo->[0]{'text'}; 
	if ( $parm{'slct_type'} eq 'blkID' ) {
		for ( my $i=1; $i<@{$alnInfo}; $i++ ) {
			defined $is_want{ $alnInfo->[$i]{'info'}[0] } or next; 
			print {$outFh} $alnInfo->[$i]{'text'}; 
		}
	} elsif ( $parm{'slct_type'} eq 'chrID' ) {
		for ( my $i=1; $i<@{$alnInfo}; $i++ ) {
			defined $is_want{ $alnInfo->[$i]{'info'}[4] } or defined $is_want{ $alnInfo->[$i]{'info'}[5] } or next; 
			print {$outFh} $alnInfo->[$i]{'text'}; 
		}
	} else {
		&stopErr("[Err] Unknown -slct_type [$parm{'slct_type'}]\n"); 
	}
	
}# mcs_blkByList() 

sub msc_glist2pairs {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	
	my ($glist_arr) = &_readInGlist($parm{'in_glist'}); # [split(/\t/, $in_line)]
	my $qr_pivot = $parm{'pivot_pat'} // '.*?' ; 
	my $qr_target = $parm{'target_pat'} // '.*?' ; 
	$qr_pivot  = qr/$qr_pivot/s; 
	$qr_target = qr/$qr_target/s; 
	
	my %used; 
	for my $ar1 (@{$glist_arr}) {
		$ar1->[0] =~ m!^Chromosome$!i and next; 
		$ar1->[0] =~ m!$qr_pivot!o or next; 
		my %list; 
		for (my $i=6; $i<@$ar1; $i++) {
			( $ar1->[$i] eq '||' or $ar1->[$i] eq '' ) and next; 
			# die "qr_target=|$qr_target|\nqr_pivot=|$qr_pivot|\n"; 
			$ar1->[$i] =~ m!$qr_target!o or next; 
			$ar1->[$i] =~ s!\-[se]$!!gi; 
			$ar1->[$i] =~ s!^[^\s:]+:!!gi; 
			$list{ $ar1->[$i] } = 1; 
		}
		$list{ $ar1->[4] } = 1; 
		my @tb = sort keys %list; 
		for (my $j=0; $j<@tb; $j++) {
			for (my $k=$j+1; $k<@tb; $k++) {
				my $tk = "$tb[$j]\t$tb[$k]"; 
				defined $used{$tk} and next; 
				$used{$tk} = 1; 
				print {$outFh} "$tk\n"; 
			}
		}
	}
}# msc_glist2pairs() 

sub msc_aln2table {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	
	my ($chr2gen, $gen2loc)  = &_readInGff($parm{'in_gff'}); 
	my ($alnInfo) = &_readInAln($parm{'in_aln'}); 
	defined $alnInfo->[0]{'info'} or shift(@{$alnInfo}); 
	
	print {$outFh} join("\t", qw/BlkID Chrom1 Start1 End1 Chrom2 Start2 End2 Strand AlnScore AlnEvalue AlnNumber Gene1 Gene2 Ka Ks KaKs AVG_Ks Med_Ks UsedKs AVG_Ka Med_Ka UsedKa AVG_KaKs Med_KaKs UsedKaKs/)."\n"; 
	for (my $i=0; $i<@{$alnInfo}; $i++) {
		my (@gen1, @gen2, @ka, @ks, @w); 
		for my $ar1 (@{$alnInfo->[$i]{'pair'}}) {
			push(@gen1, $ar1->[0]); 
			push(@gen2, $ar1->[1]); 
			push(@ka, $ar1->[3]); 
			push(@ks, $ar1->[4]); 
			push(@w, $ar1->[5]); 
		}
		my ($alnID, $score, $eval, $npair, $chr1, $chr2, $str) = @{$alnInfo->[$i]{'info'}}; 
		$str eq 'plus'  and $str = '+'; 
		$str eq 'minus' and $str = '-'; 
		my ($alnS1, $alnE1) = ( sort {$a<=>$b} ( $gen2loc->{$gen1[0]}[1], $gen2loc->{$gen1[0]}[2], $gen2loc->{$gen1[-1]}[1], $gen2loc->{$gen1[-1]}[2] ) )[0,3]; 
		my ($alnS2, $alnE2) = ( sort {$a<=>$b} ( $gen2loc->{$gen2[0]}[1], $gen2loc->{$gen2[0]}[2], $gen2loc->{$gen2[-1]}[1], $gen2loc->{$gen2[-1]}[2] ) )[0,3]; 
		print {$outFh} join("\t", $alnID, 
		 $chr1, $alnS1, $alnE1, 
		 $chr2, $alnS2, $alnE2, 
		 $str, $score, $eval, $npair, 
		 join(',', @gen1), 
		 join(',', @gen2), 
		 join(',', @ka), 
		 join(',', @ks), 
		 join(',', @w), 
		 &_avg(@ks), &_avg(@ka), &_avg(@w) 
		)."\n"; 
	}
}# msc_aln2table() 

sub msc_aln2list{
	my %parm = $ms_obj->_setHashFromArr(@_); 
	$parm{'addChr'} //= 0; 
	$parm{'srt_by'} //= 'min'; 
	
	my ($chr2gen, $gen2loc) = &_readInGff($parm{'in_gff'}); 
	my ($chr2gen_tgt, $gen2loc_tgt) = ($chr2gen, $gen2loc); 
	if (defined $parm{'tgt_gff'} and $parm{'tgt_gff'} ne '') {
		($chr2gen_tgt, $gen2loc_tgt) = &_readInGff($parm{'tgt_gff'}); 
	}
	my ($alnInfo) = &_readInAln($parm{'in_aln'}); 
	defined $alnInfo->[0]{'info'} or shift(@{$alnInfo}); 
	
	print {$outFh} join("\t", qw/Chromosome IntraDep InterDep InterTax Pivot PivotStart Blocks/)."\n"; 
	for my $chrID ( sort keys %$chr2gen ) {
		# ([genID, chrS], [], ...)
		my @glist = map { [ $_->[0], $_->[1] ] } sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @{$chr2gen->{$chrID}}; 
		# {genID} => index_number
		my %gidx; for ( my $i=0; $i<@glist; $i++ ) { $gidx{ $glist[$i][0] } = $i; } 
		my @sub_alnInfo; 
		for my $tr (@$alnInfo) {
			my ($alnID, $score, $eval, $npair, $chr1, $chr2, $str) = @{$tr->{'info'}}; 
			# When gene-1 and gene-2 locate in the same reference chrom, we only record this pair once. 
			for (my $na=0; $na<2; $na++) {
				my $nb = 1-$na; 
				if ( defined $gidx{ $tr->{'pair'}[0][$na] } and defined $gen2loc_tgt->{ $tr->{'pair'}[0][$nb] } ) {
					# The gene-$na is on the reference chrom, and gene-$nb is on the target chrom. 
					push(@sub_alnInfo, $tr); 
					my @t_sloci = map { $gidx{ $_->[$na] } } @{$tr->{'pair'}}; 
					my $ic = $ms_obj->ins_calc( \@t_sloci, 0 ); 
					$sub_alnInfo[-1]{'idx'} = $ic->{ $parm{'srt_by'} }; 
					$parm{'srt_by'} eq 'COUNT' and $sub_alnInfo[-1]{'idx'} *= -1; 
					last; 
				}
			}
		}# for my $tr (@$alnInfo) 
		@sub_alnInfo = sort { $a->{'idx'} <=> $b->{'idx'} } @sub_alnInfo; 
		my ($glist_fed, $depth, $use_tax) = &_fill_glist_wiALN(
		 'glist'=>\@glist, 
		 'alnInfo'=>\@sub_alnInfo, 
		 'gindex'=>\%gidx, 
		 'addChr'=>$parm{'addChr'}, 
		); 
		for ( my $i=0; $i<@$glist_fed; $i++ ) {
			@{$glist_fed->[$i]} = map { (defined $_) ? $_ : '' ; } @{$glist_fed->[$i]}; 
			my (@blk_tt, @blk_dd); 
			my $ref_tax = substr($chrID, 0, 2); 
			my $intra_dep = 0; 
			if ( defined $depth->[$i] ) {
				for my $blk_tax ( sort keys %{$depth->[$i]} ) {
					if ( $blk_tax eq $ref_tax ) {
						$intra_dep = $depth->[$i]{$blk_tax}; 
					} else {
						push(@blk_dd, $depth->[$i]{$blk_tax}); 
						push(@blk_tt, $blk_tax); 
					}
				}
			}
			scalar(@blk_dd) == 0 and do { @blk_dd=(0); @blk_tt=(''); }; 
			print {$outFh} join("\t", $chrID, $intra_dep, join(',', @blk_dd), join(',', @blk_tt), @{$glist_fed->[$i]})."\n"; 
		}
	}# for my $chrID ( sort keys %$chr2gen ) 
}# msc_aln2list() 

####################################################################################################
# Internal sub-functions that should only be used by other sub-functions. 
####################################################################################################

# Input  : &_fill_glist_wiALN( 'glist'=>\@glist, 'alnInfo'=>\@sub_alnInfo, 'gindex'=>\%gidx, 'addChr'=>$parm{'addChr'} ); 
# Return : (\@glist_wiALN, \@depth, \%use_tax)
#  @glist_wiALN = ( [genID, chrS, ] )
sub _fill_glist_wiALN {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	$parm{'addChr'} //= 0; 
	$parm{'srt_by'} //= 'min'; 
	my @glist_wiALN = @{$parm{'glist'}}; 
	unless ( defined $parm{'gindex'} ) {
		for ( my $i=0; $i<@glist_wiALN; $i++ ) { $parm{'gindex'}{ $glist_wiALN[$i][0] } = $i; }
	}
	my (@depth, %use_tax); 
	
	ALN: for my $tr ( @{$parm{'alnInfo'}} ) {
		my ($alnID, $score, $eval, $npair, $chr1, $chr2, $str) = @{$tr->{'info'}}; 
		defined $chr1 or &stopErr("[Err] chr1 not defined in [@{$tr->{'info'}}]\n"); 
		my @chr = ($chr1, $chr2); 
		my @gpair = @{$tr->{'pair'}}; 
		for ( my $na=0; $na<2; $na++ ) {
			defined $parm{'gindex'}{$gpair[0][$na]} or next; 
			my $nb = 1-$na; 
			my @tax = ( substr($chr1, 0, 2), substr($chr2, 0, 2) ); 
			my ( $minIdx, $maxIdx ) = ( $parm{'gindex'}{ $gpair[0][$na] }, $parm{'gindex'}{ $gpair[-1][$na] } ); 
			$minIdx > $maxIdx and do { ($minIdx, $maxIdx) = ($maxIdx, $minIdx); @gpair = reverse(@gpair);  }; 
			# Find the smallest track No. (column_index) for current alnID; 
			my $track = 1; # The 0-index column stores 'genID' information. 
			for ( ;1;$track++ ) {
				my $is_ok = 1; 
				for (my $j=$minIdx; $j<=$maxIdx; $j++) {
					defined $glist_wiALN[$j][$track] and $glist_wiALN[$j][$track] ne '' and do { $is_ok=0; last; }; 
				}
				$is_ok == 1 and last; 
			}
			for ( my ($j,$k)=($minIdx, 0); $j<=$maxIdx and $k<@gpair; $j++ ) {
				$depth[$j]{$tax[$nb]} ++; 
				defined $use_tax{$tax[$nb]} or $use_tax{$tax[$nb]} = 1; 
				$glist_wiALN[$j][0] ne $gpair[$k][$na] and do { $glist_wiALN[$j][$track] = '||'; next; }; 
				$glist_wiALN[$j][$track] = $gpair[$k][$nb]; 
				$k == 0       and $glist_wiALN[$j][$track] .= "-s"; 
				$k == $#gpair and $glist_wiALN[$j][$track] .= "-e"; 
				$parm{'addChr'} and $glist_wiALN[$j][$track] = "$chr[$nb]:$glist_wiALN[$j][$track]"; 
				$k++; 
			}
		}
	}# ALN: for my $tr ( @{$parm{'alnInfo'}} ) 
	return (\@glist_wiALN, \@depth, \%use_tax); 
}# _fill_glist_wiALN()

sub _readInGlist {
	my $inFh = &openFH(shift, '<'); 
	my @back_arr; 
	while (<$inFh>) {
		chomp; 
		push(@back_arr, [split(/\t/, $_)]); 
	}
	return \@back_arr; 
}# _readInGlist() 

# Return : (\%chr2gen, \%gen2loc)
#  %chr2gen : {chromID} => [ [genID, chrS, chrE], [], ... ]
#  %gen2loc : {genID} => [chrID, chrS, chrE]
sub _readInGff {
	my $inFh = &openFH(shift, '<'); 
	my (%chr2gen, %gen2loc); 
	while (<$inFh>) {
		chomp; 
		my ( $chrID, $genID, $chrS, $chrE ) = split(/\t/, $_); 
		push( @{$chr2gen{$chrID}}, [$genID, $chrS, $chrE] ); 
		defined $gen2loc{$genID} and &tsmsg("[Msg] genID [$genID] repeats, and loc infor masked by latter one\n"); 
		$gen2loc{$genID} = [ $chrID, $chrS, $chrE ]; 
	}
	close($inFh); 
	return (\%chr2gen, \%gen2loc); 
}# _readInGff() 
# Return : (\@alnInfo)
#  @alnInfo : 
#    [idx_num]{'info'} => [alnID, score, evalue, Num_pairs, chrID_1, chrID_2, strand(plus|minus)]
#    [idx_num]{'pair'} => [ [genID_1, genID_2, pair_evalue, Ka, Ks, Ka/Ks], [], ... ]
#    [idx_num]{'text'} => $text_of_current_block
sub _readInAln{
	my $inFh = &openFH(shift, '<'); 
	my @alnInfo; 
	my $aln_idx = 0; 
	while (<$inFh>) {
		chomp; 
		if ( m!^\s*$|^####|^# \S+:\s+\S+$|# (MAX GAPS|Number of collinear genes|Number of all genes):! ) {
			$alnInfo[$aln_idx]{'text'} .= "$_\n"; 
			next; 
		}
		if ( m!^## Alignment ! ) {
			# ## Alignment 0: score=7137.0 e_value=0 N=147 Ma1&Ma1 plus
			m!## Alignment (\d+): score=(\S+) e_value=(\S+) N=(\d+) ([^\&\s]+)\&([^\&\s]+) (plus|minus|X+)! or die "ALN:$_\n"; 
			my ($alnID, $score, $eval, $n, $chr1, $chr2, $str)
			=  ($1,     $2,     $3,    $4, $5,    $6,    $7); 
			$aln_idx++; 
			$alnInfo[$aln_idx]{'info'} = [$alnID, $score, $eval, $n, $chr1, $chr2, $str]; 
			$alnInfo[$aln_idx]{'text'} .= "$_\n"; 
		} elsif ( m!^\s*(\d+)\-\s*(\d+):\t(\S+)\t(\S+)\s*(\S+)(?:\t(\S+)\t(\S+)(?:\t(\S+))?)?$! ) { 
			# #  0-  0:        Cma_000007      Cma_000973        2e-57
			$aln_idx > -1 or &stopErr("Too early to line: $_\n"); 
			my ($alnID, $alnID_id, $gid1, $gid2, $eval, $tka, $tks, $tw) = ($1, $2, $3, $4, $5, $6, $7, $8); 
			$tka //= ''; $tks //= ''; $tw //= ''; 
			$alnID == $alnInfo[$aln_idx]{'info'}[0] or &stopErr("[Err] line_alnID=$alnID not fitting upper level (alnID=$alnInfo[$aln_idx]{'info'}[0]).\n"); 
			push(@{$alnInfo[$aln_idx]{'pair'}}, [$gid1, $gid2, $eval, $tka, $tks, $tw]); 
			$alnInfo[$aln_idx]{'text'} .= "$_\n"; 
		} else {
			&stopErr("[Err] Unable to parse line: $_\n"); 
		}
	}
	close($inFh); 
	return(\@alnInfo); 
}# _readInAln

# Input  : .bn6 file (blastn -outfmt 6)
# Return : (\@bn6Info)
#   [idx]{'arr'} = [split(/\t/, $_)]
#   [idx]{'k2v'}{$key} = $value; 
#     $key could be qw/qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand/
sub _readInBn6 {
	my $inFh = &openFH(shift, '<'); 
	my @bn6Info; 
	while (<$inFh>) {
		chomp; m!^\s*$! and next; 
		my @ta = split(/\t/, $_); 
		push(@bn6Info, {}); 
		$bn6Info[-1]{'arr'} = \@ta; 
		$bn6Info[-1]{'k2v'}{'qseqid'}    = $ta[0]; 
		$bn6Info[-1]{'k2v'}{'sseqid'}    = $ta[1]; 
		$bn6Info[-1]{'k2v'}{'pident'}    = $ta[2]; 
		$bn6Info[-1]{'k2v'}{'length'}    = $ta[3]; 
		$bn6Info[-1]{'k2v'}{'mismatch'}  = $ta[4]; 
		$bn6Info[-1]{'k2v'}{'gapopen'}   = $ta[5]; 
		$bn6Info[-1]{'k2v'}{'qstart'}    = $ta[6]; 
		$bn6Info[-1]{'k2v'}{'qend'}      = $ta[7]; 
		$bn6Info[-1]{'k2v'}{'sstart'}    = $ta[8]; 
		$bn6Info[-1]{'k2v'}{'send'}      = $ta[9]; 
		$bn6Info[-1]{'k2v'}{'evalue'}    = $ta[10]; 
		$bn6Info[-1]{'k2v'}{'bitscore'}  = $ta[11]; 
		$bn6Info[-1]{'k2v'}{'qlen'}      = $ta[12] // ''; 
		$bn6Info[-1]{'k2v'}{'slen'}      = $ta[13] // ''; 
		$bn6Info[-1]{'k2v'}{'sstrand'}   = $ta[14] // ''; 
		$bn6Info[-1]{'k2v'}{'bitscore'} =~ s!^\s+|\s+$!!g; 
	}
	return \@bn6Info; 
}# _readInBn6() 

sub _avg {
	my ($avg, $med) = ('nan', 'nan'); 
	my @unn = grep { defined $_ and $_ !~ m!^(|\-?nan|\-?inf)$!i and $_ >= 0 } @_; 
	my $validNum = scalar(@unn); 
	if ( $validNum > 0 ) {
		my $ic = $ms_obj->ins_calc(\@unn, 0); 
		$avg = $ic->{'MEAN'}; 
		$med = $ic->{'MEDIAN'}; 
	}
	return ($avg, $med, $validNum); 
}# _avg() 

sub _write_td {
	my $color = ''; 
	defined $_[1] and $color = " bgcolor='$_[1]'"; 
	return "<td$color>$_[0]</td>"; 
}# _write_td() 

# Input  : _index_glist( $chr2gen )
# Return : \%gindex
#  {$genID} => increasing_1by1_number; 
sub _index_glist {
	my $chr2gen = shift; 
	my %gidx; 
	my $ii = 0; 
	for my $chrID (sort keys %$chr2gen) {
		for my $ar1 (sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @{$chr2gen->{$chrID}}) {
			defined $gidx{$ar1->[0]} and &stopErr("[Err] Repeat genID [$ar1->[0]]\n"); 
			$gidx{$ar1->[0]} = $ii; 
			$ii++; 
		}
	}
	return \%gidx; 
}# _index_glist() 

# Input  : _filter_eval( $bn6Info, $max_eval )
# Return : \@filtered_bn6Info 
sub _filter_eval {
	my ($bn6Info, $max_eval) = @_; 
	( defined $max_eval and $max_eval >= 0) or return $bn6Info; 
	my @back_bn6Info; 
	for my $ar1 (@$bn6Info) {
		$ar1->{'k2v'}{'evalue'} <= $max_eval or next; 
		push(@back_bn6Info, $ar1); 
	}
	return \@back_bn6Info; 
}# _filter_eval() 

# Input  : _tandemGrp( 'id_grp'=>$id_grp, 'max_distN'=>$parm{'max_tandN'}, 'gen2loc'=>$gen2loc, 'gindex'=>\%gidx,  )
# Return : \@tand_gene 
#   ([blk0_gen1, blk0_gen2, ...], ...)
sub _tandemGrp {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	$parm{'max_distN'} //= 2; 
	unless ( defined $parm{'gindex'} and keys %{$parm{'gindex'}} > 0 ) {
		my $i=0; 
		for my $genID ( sort { $parm{'gen2loc'}->{$a}[0] cmp $parm{'gen2loc'}->{$b}[0] 
		             || $parm{'gen2loc'}{$a}[1] <=> $parm{'gen2loc'}{$b}[1] 
					 || $parm{'gen2loc'}{$a}[2] <=> $parm{'gen2loc'}{$b}[2] 
		      } keys %{$parm{'gen2loc'}} 
		) {
			$parm{'gindex'}{$genID} = $i; 
			$i ++; 
		}
	}
	my %tand_relate; 
	# my @tand_gene; 
	for my $ar1 (@{$parm{'id_grp'}}) {
		my @ta = sort { $parm{'gindex'}{$a} <=> $parm{'gindex'}{$b} } @$ar1; 
		for (my $i=0; $i+1<@ta; $i++) {
			my ($g1, $g2) = (@ta[$i, $i+1]); 
			$parm{'gen2loc'}{$g1}[0] eq $parm{'gen2loc'}{$g2}[0] or next; 
			$parm{'gindex'}{$g2} - $parm{'gindex'}{$g1} <= $parm{'max_distN'} or next; 
			$tand_relate{$g1}{$g2} = 1; 
			$tand_relate{$g2}{$g1} = 1; 
			# if ( scalar(@tand_gene) > 0 ) {
			# 	if ( $tand_gene[-1][-1] eq $g1 ) {
			# 		push(@{$tand_gene[-1]}, $g2); 
			# 	} else {
			# 		push(@tand_gene, [$g1, $g2]); 
			# 	}
			# } else {
			# 	@tand_gene = ([$g1, $g2]); 
			# }
		}
	}
	my $tand_gene = &_chain_grp(\%tand_relate); 
	return $tand_gene; 
}# _tandemGrp() 

# Input  : _filter_cscore( $bn6Info, $min_cscore, $gen2loc ) 
# Return : \@filtered_bn6Info
sub _filter_cscore {
	my ( $bn6Info, $min_cscore, $gen2loc ) = @_; 
	( defined $min_cscore and $min_cscore > 0 ) or return $bn6Info; 
	my @back_bn6Info; 
	if ( defined $gen2loc ) {
		my %gen2tax; 
		my %taxpair2idx; 
		for my $genID ( keys %$gen2loc ) {
			my $taxID = substr( $gen2loc->{$genID}[0], 0, 2 ); 
			$gen2tax{ $genID } = $taxID; 
		}
		for (my $i=0; $i<@$bn6Info; $i++ ) {
			my $ar1 = $bn6Info->[$i]; 
			my $tax1 = $gen2tax{ $ar1->{'k2v'}{'qseqid'} } // 'NA'; 
			my $tax2 = $gen2tax{ $ar1->{'k2v'}{'sseqid'} } // 'NA'; 
			my $taxpair = join( "\t", sort ($tax1, $tax2) ); 
			push(@{$taxpair2idx{$taxpair}}, $i); 
		}
		my @back_idx; 
		for my $taxpair ( keys %taxpair2idx ) {
			my @tt_bn6Info = @{$bn6Info}[ @{$taxpair2idx{$taxpair}} ]; 
			my $best_score = &_bestScore( \@tt_bn6Info ); 
			for my $ii ( @{$taxpair2idx{$taxpair}} ) {
				my $ar1 = $bn6Info->[$ii]; 
				my $cscore = $ar1->{'k2v'}{'bitscore'} / &mathSunhh::max( $best_score->{ $ar1->{'k2v'}{'qseqid'} }, $best_score->{ $ar1->{'k2v'}{'sseqid'} } ); 
				$cscore >= $min_cscore or next; 
				push(@back_idx, $ii); 
			}
		}
		for my $ii ( sort { $a<=>$b } @back_idx ) {
			push( @back_bn6Info, $bn6Info->[$ii] ); 
		}
	} else {
		my $best_score = &_bestScore( $bn6Info ); 
		for my $ar1 ( @$bn6Info ) {
			my $cscore = $ar1->{'k2v'}{'bitscore'} / &mathSunhh::max( $best_score->{ $ar1->{'k2v'}{'qseqid'} }, $best_score->{ $ar1->{'k2v'}{'sseqid'} } ); 
			$cscore >= $min_cscore or next; 
			push(@back_bn6Info, $ar1); 
		}
	}
	
	return \@back_bn6Info; 
}# _filter_cscore() 

# Input  : _id2net(\%id2id, $seedID, \%inNetID)
# Function : Retrieve all IDs that can be net linked to $seedID by \%id2id; 
# Return : \%inNetID
sub _id2net {
	my ($id2id, $seedID, $inNetID) = @_; 
	$inNetID //= {}; 
	$inNetID->{ $seedID } = 1; 
	my @nextID = grep { !defined $inNetID->{$_} } keys %{$id2id->{$seedID}}; 
	if ( scalar(@nextID) > 0 ) {
		for my $nID (@nextID) {
			defined $inNetID->{$nID} and next; 
			&_id2net($id2id, $nID, $inNetID); 
		}
	} else {
		return $inNetID; 
	}
	return $inNetID; 
}# _id2net

# Input  : _group_id( $bn6Info ) 
# Return : \@id_grp
#  [[grp0_gen1, grp0_gen2, ...], [grp1_gen1, grp2_gen2, ...], ...]
sub _group_id {
	my ($bn6Info) = @_; 
	my %id2id; 
	for my $ar1 ( @$bn6Info ) {
		my $g1 = $ar1->{'k2v'}{'qseqid'}; 
		my $g2 = $ar1->{'k2v'}{'sseqid'}; 
		$id2id{$g1}{$g2} ++; 
		$id2id{$g2}{$g1} ++; 
	}
	# my ($grps) = &_chain_grp(\%id2id); 
	my ($grps) = &_1step_grp(\%id2id); 
	return $grps; 
}# _group_id() 

sub _1step_grp {
	my $id2id = shift; 
	my @grps; 
	for my $id1 (sort keys %$id2id) {
		push(@grps, [sort keys %{$id2id->{$id1}}]); 
	}
	return \@grps; 
}# _1step_grp() 

sub _chain_grp{
	my $id2id = shift; 
	my %inGrp; 
	my @grps; 
	for my $id (keys %$id2id) {
		defined $inGrp{$id} and next; 
		my ($all_net2id) = &_id2net($id2id, $id); 
		push(@grps, [ sort keys %$all_net2id ]); 
		for my $id2 ( @{$grps[-1]} ) {
			$inGrp{$id2} = 1; 
		}
	}
	return \@grps; 
}# _chain_grp() 

# Input  : _bestScore( $bn6Info ) 
# Return : \%best_score
#  {$genID} => max_score_within_q_and_s_paired_blast
sub _bestScore {
	my ( $bn6Info ) = @_; 
	my %tax; 
	my %best_score; 
	for my $ar1 (@$bn6Info) {
		$best_score{ $ar1->{'k2v'}{'qseqid'} } //= $ar1->{'k2v'}{'bitscore'}; 
		$best_score{ $ar1->{'k2v'}{'sseqid'} } //= $ar1->{'k2v'}{'bitscore'}; 
		$best_score{ $ar1->{'k2v'}{'qseqid'} } < $ar1->{'k2v'}{'bitscore'} and $best_score{ $ar1->{'k2v'}{'qseqid'} } = $ar1->{'k2v'}{'bitscore'}; 
		$best_score{ $ar1->{'k2v'}{'sseqid'} } < $ar1->{'k2v'}{'bitscore'} and $best_score{ $ar1->{'k2v'}{'sseqid'} } = $ar1->{'k2v'}{'bitscore'}; 
	}
	return \%best_score; 
}# _bestScore() 

sub _filter_tand {
	my ($bn6Info, $max_tandN, $gen2loc) = @_; 
	my ($id_grp) = &_group_id( $bn6Info ); 
	my $tand_gene = &_tandemGrp(
	 'id_grp' => $id_grp, 
	 'max_distN' => $max_tandN, 
	 'gen2loc' => $gen2loc, 
	); 
	my %is_tand; 
	for my $ar1 ( @$tand_gene ) {
		map { $is_tand{ $_ } = 1; } @$ar1; 
	}
	my @back_bn6Info; 
	for my $ar1 (@$bn6Info) {
		$is_tand{ $ar1->{'arr'}[0] } and next; 
		$is_tand{ $ar1->{'arr'}[1] } and next; 
		push(@back_bn6Info, $ar1); 
	}
	return \@back_bn6Info; 
}# _filter_tand() 

sub _filter_repeat {
	my ($bn6Info, $repeat_eval) = @_; 
	$repeat_eval //= 0.05; 
	
	my %pair; 
	my $cnt_pair = 0; 
	for my $ar1 (@$bn6Info) {
		my ($g1, $g2) = ($ar1->{'arr'}[0], $ar1->{'arr'}[1]); 
		defined $pair{$g1}{$g2} or defined $pair{$g2}{$g1} or $cnt_pair ++; 
		$pair{$g1}{$g2} = 1; 
		$pair{$g2}{$g1} = 1; 
	}
	my $cnt_genes = scalar(keys %pair); 
	my %cnt_g; 
	for (keys %pair) {
		$cnt_g{$_} = scalar(keys %{$pair{$_}}); 
	}
	my @back_bn6Info; 
	for my $ar1 ( @$bn6Info ) {
		my ($g1, $g2) = ($ar1->{'arr'}[0], $ar1->{'arr'}[1]); 
		my $adj_eval = $ar1->{'k2v'}{'evalue'} ** ( ($cnt_pair/$cnt_genes) / ($cnt_g{$g1}+$cnt_g{$g2}) ); 
		$adj_eval > $repeat_eval and next; 
		push(@back_bn6Info, $ar1); 
	}
	return \@back_bn6Info; 
}# _filter_repeat() 

