#!/usr/bin/perl -w
# 20200806: For -add_KaKs function, fix the problem when the input CDS IDs have too many characters (>= 32). 
# 20200820: Fix -add_KaKs when there are duplicated gene names. 
# [20220310] There is a bug in add_KaKs function which causes no KaKs attached to the original .collinearity (-in_aln) file when using "-ncpu 1" instead of multiple cpus. I'll fix it later when I have time.
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 
use mcsSunhh; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 

use Parallel::ForkManager;


use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
 "in_gff:s", # in_gff is similar to mcscan's .bed file, with format as : (OO)N+ \t GeneID \t Start \t End 
 "in_aln:s", # This is the output ".collinearity" file of mcscanX, which should be similar to .align of mcscan or previous mcscanX version. 
 "in_blast:s", # This is the input ".blast" file of mcscanX, 
 "in_table:s", # This is the output of -aln2table. 
 "out:s", # output filename 
 "aln2list!", "addChr!", "tgt_gff:s", "srt_by:s", "raw_order!", 
 "aln2table!",
   "useYN!", 
 "table2aln!", # require '-in_table'
 "sepTable!", 
   "sepColN:s", 
 "aln2pairList!", 
 "glist2pairs!", "in_glist:s", "pivot_pat:s", "target_pat:s", 
 "slctBlk!", "slct_list:s", "slct_colN:i", "slct_type:s", 
 "glist2html!", "pivot_chrID:s", 
 "classDupGene!", "max_eval:f", "min_cscore:f", "max_proxN:i", "max_tandN:i", "tand_list:s", "only_direct!", "prox_list:s", 
 "filterBlast!", "rm_repeat!", "repeat_eval:f", "rm_tandem!", 
 "add_KaKs!", "in_pair_list:s", "fas_cds:s", "fas_prot:s", "out_pref:s","alnMethod:s", "ncpu:i", "kaks_tab:s", 
 "cvt_ctg2scf!", "scf_agp:s", 
 "trim_block!", "trim_blkSE:s", 
 "filter_blk!", "min_blkPair:i", 
 "exe_python:s", # python
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
#   -tgt_gff      [target.gff] Provide source of blocks. Only blocks with genes from this gff aligned to the backbone genes are kept. 
#   -srt_by       ['min'] Sort the gene blocks by 'min|Q1|Q3|interval_mean|COUNT'
#   -raw_order    [Boolean] Do not sort the input gff according to their locations. 
#                   Sometimes the orders of genes are different in scaffold and chromosome if there is a gene included by another one. (Rarely happen)
# 
# -aln2pairList   [Boolean] get gene pair list from aligned blocks. 
#                   Need -in_aln 
#                   May use -useYN 
# -aln2table      [Boolean] Reformat aligned blocks into one line. 
#                   Need -in_gff , -in_aln 
#                   Headers: BlkID / Chrom1 / Start1 / End1 / Chrom2 / Start2 / End2 / Strand / AlnScore / AlnEvalue / AlnNumber / Gene1 / Gene2 / 
#                            Ka / Ks / KaKs / AVG_Ks / Med_Ks / UsedKs / AVG_Ka / Med_Ka / UsedKa / AVG_KaKs / Med_KaKs / UsedKaKs
#   -useYN        [Boolean] By default, Nei-Gojobori (NG) will be chose in block Ks calculation. But
#                   I will switch to use Yang-Nielson (YN) results if given this parameter. 
#
# -sepTable       [Boolean] Separate table into one gene pair per line format. This is useful when checking gene pair's kaks. 
#                   Need -in_table 
#   -sepColN      [0,11,12,13,14,15] Only use these columns; 
# 
# -table2aln      [Boolean] requires -in_table
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
#   -prox_list    [outFileName] 
#   -only_direct  [Boolean] Only count homologous relationship between directly aligned gene pairs in .blast file. 
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
#   -kaks_tab     [NULL] Do not calculate KaKs again if given. 
#
# -cvt_ctg2scf    [Boolean] Convert contig IDs into scaffold IDs according to AGP file. 
#                   Need -in_aln ; 
#                   Please note that a contig should be related to only ONE scaffold here!!! 
#   -scf_agp      [filename] AGP file from contigs to scaffolds. 
#
# -trim_block     [Boolean] Trim the blocks according to given blocks' start-end (SE) reference genes. 
#                   Need -in_aln ; 
#                   Please note that the e_value will not be changed. 
#   -trim_blkSE   [filename] A file telling the new start and end genes in a block. 
#                   Format : 'Block_ID \\t Start_Gene \\t End_Gene'
#                   Blocks not mentioned will be output intactly. 
#                   If the format is : 'Block_ID \\t "Gene1_prev : Gene1 ; Gene2 ; Evalue ..." \\t .Add'
#                   I will add the 'Gene1-Gene2' pair to Block_ID directly after Gene1_prev. 
#                   If the Gene1_prev == '.Top', these pairs will be added to the top of the block. 
#                   If the format is : 'Block_ID \\t Start_Gene \\t .Remove'
#                   I will remove the 'Start_Gene' in related block. 
#                     The order of modification is '.Add' => '.Remove' => trimming. 
#   -useYN        [Boolean] By default, Nei-Gojobori (NG) will be chose in block Ks calculation. But
#                   I will switch to use Yang-Nielson (YN) results if given this parameter. 
#
# -filter_blk   [Boolean] 
#   -min_blkPair  [Number] Default 2. At least this number of pairs in this block. 
#
# -exe_python   [python] Could try '/usr/bin/python'
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
$opts{'useYN'}  //= 0; 
$opts{'only_direct'} //= 0; 
$opts{'exe_python'} //= 'python'; 


####################################################################################################
# Call different sub-funcitons. 
####################################################################################################

if ( $opts{'aln2list'} ) {
	&msc_aln2list(
	 'in_gff'=>$opts{'in_gff'}, 
	 'in_aln'=>$opts{'in_aln'}, 
	 'addChr'=>$opts{'addChr'}, 
	 'tgt_gff'=>$opts{'tgt_gff'}, 
	 'srt_by'=>$opts{'srt_by'}, 
	 'raw_order'=>$opts{'raw_order'}, 
	); 
} elsif ( $opts{'aln2table'} ) {
	&msc_aln2table(
	 'in_gff'=>$opts{'in_gff'}, 
	 'in_aln'=>$opts{'in_aln'}, 
	); 
} elsif ( $opts{'table2aln'} ) {
  &mcs_table2aln(
   'in_table' => $opts{'in_table'}
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
	 'prox_list'=>$opts{'prox_list'}, 
	 'only_direct' => $opts{'only_direct'}, 
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
	 'kaks_tab' => $opts{'kaks_tab'}, 
	); 
} elsif ( $opts{'cvt_ctg2scf'} ) {
	&cvt_ctg2scf_byAGP(
	 'in_aln'  => $opts{'in_aln'}, 
	 'scf_agp' => $opts{'scf_agp'}, 
	); 
} elsif ( $opts{'trim_block'} ) {
	&trim_block(
	 'in_aln'  => $opts{'in_aln'},
	 'trim_blkSE' => $opts{'trim_blkSE'}, 
	);  
} elsif ( $opts{'filter_blk'} ) {
	&filter_block(
	 'in_aln'  => $opts{'in_aln'}, 
	 'min_blkPair' => $opts{'min_blkPair'}, 
	); 
} elsif ( $opts{'aln2pairList'} ) {
	&mcs_aln2pair (
	 'in_aln' => $opts{'in_aln'}, 
	 'useYN' => $opts{'useYN'}, 
	); 
} elsif ( $opts{'sepTable'} ) {
	$opts{'sepColN'} //= '0,11,12,13,14,15'; 
	&sep_tab (
	  'in_table'=> $opts{'in_table'}, 
	  'sepColN' => $opts{'sepColN'}, 
	); 
} else {
	&usage(); 
}
####################################################################################################
# Sub-functions to be called 
####################################################################################################

sub sep_tab {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	$parm{'sepColN'} //= '0,11,12,13,14,15'; 
	my @CV = split(/,/, $parm{'sepColN'}); 
	open F,'<',"$parm{'in_table'}" or die; 
	while (<F>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		if ( $ta[0] eq 'BlkID' ) {
			print STDOUT join("\t", @ta[@CV])."\n"; 
			next; 
		}
		my $maxN = -1; 
		for my $tv (@CV) {
			my @tc = split(/,/, $ta[$tv]);
			$maxN < scalar(@tc) and $maxN = scalar(@tc);
		}
		for (my $i=0; $i<$maxN; $i++) {
			my @o; 
			for (my $j=0; $j<@CV; $j++) {
				my $tv = $CV[$j];
				my @tc = split(/,/, $ta[$tv]);
				if ( defined $tc[$i] ) {
					push(@o, $tc[$i]);
				} else {
					push(@o, $tc[0]);
				}
			}
			print STDOUT join("\t", @o)."\n"; 
		}
	}
	close F; 
}# sep_tab () 

sub mcs_addKaKs {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	$parm{'fas_prot'} //= ''; 
	$parm{'out_pref'} //= 'paired'; 
	$parm{'ncpu'} //= 1; 
	if ( defined $parm{'kaks_tab'} and $parm{'kaks_tab'} ne '' ) {
		&attach_KaKs( $parm{'in_aln'}, $parm{'kaks_tab'}, "$parm{'in_aln'}.ks" ); 
		return; 
	}
	
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
		my ($alnInfo) = &_readInAln( $parm{'in_aln'}, $opts{'useYN'} ); 
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
		&attach_KaKs( $parm{'in_aln'}, "$parm{'out_pref'}.ks", "$parm{'in_aln'}.ks", $header, $pair2ks ); 
	}
	return ; 
}# mcs_addKaKs() 


sub attach_KaKs {
	my ($alnF, $kaksF, $outF, $header, $pair2ks) = @_; 
	my ($alnInfo) = &_readInAln( $alnF, $opts{'useYN'} );
	if (!defined $header or !defined $pair2ks) {
		($header, $pair2ks) = &_readInKsTab( $kaksF );
	}
	my $oksFh = &openFH( $outF, '>' ); 
	for my $ar1 (@$alnInfo) {
		for my $ln (split(/\n/, $ar1->{'text'})) {
			chomp($ln); 
			if ( $ln =~ m/^\s*(#|$)/ ) {
				print {$oksFh} "$ln\n"; 
				next; 
			}
			my @ta = split(/\t/, $ln); 
			$pair2ks->{$ta[1]}{$ta[2]} //= [ 'nan', 'nan', 'nan', 'nan' ]; 
			print {$oksFh} join( "\t", @ta[0..3], @{ $pair2ks->{$ta[1]}{$ta[2]} } )."\n"; 
		}
	}
	close ($oksFh); 
	return; 
}# sub attach_KaKs() 

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
	# Edit on 2020-08-06 to change input sequence ID to a/b before running jcvi; 
	my %parm = $ms_obj->_setHashFromArr(@_); 
	$parm{'fas_prot'} //= ''; 
	$parm{'out_pref'} //= 'paired'; 
	$parm{'alnMethod'} //= 'muscle'; 
	defined $parm{'fas_cds'} or &stopErr("[Err] -fas_cds must be given.\n"); 
	defined $parm{'in_pair_list'} or &stopErr("[Err] -in_pair_list is needed.\n"); 
	my %cds_seq = %{ $fs_obj->save_seq_to_hash('faFile'=>$parm{'fas_cds'}) }; 
	my %pep_seq; 
	my $ofh_cds = &openFH("$parm{'out_pref'}.cds.fasta", '>'); 
	my $ofh_pep; 
	my $med_prot = ''; 
	if ( $parm{'fas_prot'} ne '' ) {
		$med_prot = "$parm{'out_pref'}.cds.fasta.pep"; 
		$ofh_pep = &openFH($med_prot, '>'); 
		%pep_seq = %{ $fs_obj->save_seq_to_hash('faFile'=>$parm{'fas_prot'}) }; 
	}
	my (%id_old2new, $tmpCnt1, %id_new2old); 
	for my $tmpLID (map { @{$_}[0,1] } &fileSunhh::load_tabFile($parm{'in_pair_list'})) {
		defined $cds_seq{$tmpLID} or &stopErr("[Err] Failed to find sequence for ID [$tmpLID]\n"); 
		#if ( defined $id_old2new{$tmpLID} ) {
		#	$tmpCnt1 = $id_old2new{$tmpLID}; 
		#} else {
		#	$tmpCnt1 ++; 
		#	$id_old2new{$tmpLID} = $tmpCnt1; 
		#	$id_new2old{$tmpCnt1} = $tmpLID; 
		#}
		$tmpCnt1 ++; 
		push(@{$id_old2new{$tmpLID}}, $tmpCnt1); 
		$id_new2old{$tmpCnt1} = $tmpLID; 
		$cds_seq{$tmpLID}{'seq'} =~ s!\s!!g; 
		$cds_seq{$tmpLID}{'seq'} =~ s!(.{60})!$1\n!g; chomp($cds_seq{$tmpLID}{'seq'}); 
		print {$ofh_cds} ">$tmpCnt1\n$cds_seq{$tmpLID}{'seq'}\n"; 
		if ($med_prot ne '') {
			defined $pep_seq{$tmpLID} or &stopErr("[Err] Failed to find pep sequence for ID [$tmpLID]\n"); 
			$pep_seq{$tmpLID}{'seq'} =~ s!\s!!g; 
			$pep_seq{$tmpLID}{'seq'} =~ s!(.{60})!$1\n!g; chomp($pep_seq{$tmpLID}{'seq'}); 
			print {$ofh_pep} ">$tmpCnt1\n$pep_seq{$tmpLID}{'seq'}\n"; 
		}
	}
	close($ofh_cds); 
	close($ofh_pep); 
	# &exeCmd_1cmd("$opts{'exe_python'} -m jcvi.apps.ks prepare $parm{'in_pair_list'} $parm{'fas_cds'} $parm{'fas_prot'} -o $parm{'out_pref'}.cds.fasta") and &stopErr("[Err] Failed prepare\n"); 
	&exeCmd_1cmd("$opts{'exe_python'} -m jcvi.apps.ks calc    --msa=$parm{'alnMethod'} $med_prot $parm{'out_pref'}.cds.fasta -o $parm{'out_pref'}.cds.fasta.ks"); 
	my ( $header, $pair2ks ) = &_readInKsLis( "$parm{'out_pref'}.cds.fasta.ks" ); 
	my $paFh = &openFH($parm{'in_pair_list'}, '<'); 
	my $oksFh = &openFH("$parm{'out_pref'}.ks", '>'); 
	# print {$oksFh} join("\t", qw/gene1 gene2/, @$header[1 .. $#{$header}])."\n"; 
	print {$oksFh} join("\t", @$header[0 .. $#{$header}])."\n"; 
	while (<$paFh>) {
		m/^\s*(#|$)/ and next; 
		chomp; s!\s+$!!; 
		my @ta = split(/\t/, $_); 
		my ($g1Ori, $g2Ori) = @ta[0,1]; 
		my $find_ks = 0; 
		FFF1:
		for my $g1 (@{$id_old2new{$g1Ori}}) {
			for my $g2 (@{$id_old2new{$g2Ori}}) {
				defined $pair2ks->{$g1}{$g2} or next; 
				$find_ks = 1; 
				print {$oksFh} join("\t", $g1Ori, $g2Ori, @{$pair2ks->{$g1}{$g2}})."\n"; 
				last FFF1; 
			}
		}
		$find_ks == 0 and do { &tsmsg("[Err] No Ks result found for line: $_\n"); next; }; 
	}
	close ($paFh); 
	
	return ($header, $pair2ks);
}# _lis2ks() 

sub _lis2ks_ori {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	$parm{'fas_prot'} //= ''; 
	$parm{'out_pref'} //= 'paired'; 
	$parm{'alnMethod'} //= 'muscle'; 
	defined $parm{'fas_cds'} or &stopErr("[Err] -fas_cds must be given.\n"); 
	defined $parm{'in_pair_list'} or &stopErr("[Err] -in_pair_list is needed.\n"); 
	
	&exeCmd_1cmd("$opts{'exe_python'} -m jcvi.apps.ks prepare $parm{'in_pair_list'} $parm{'fas_cds'} $parm{'fas_prot'} -o $parm{'out_pref'}.cds.fasta") and &stopErr("[Err] Failed prepare\n"); 
	my $med_prot = ( $parm{'fas_prot'} eq '' ) ? "" : "$parm{'out_pref'}.cds.fasta.pep" ; 
	&exeCmd_1cmd("$opts{'exe_python'} -m jcvi.apps.ks calc    --msa=$parm{'alnMethod'} $med_prot $parm{'out_pref'}.cds.fasta -o $parm{'out_pref'}.cds.fasta.ks"); 
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
}# _lis2ks_ori() 


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
	$parm{'only_direct'} //= 0; 
	
	my ($chr2gen, $gen2loc) = &_readInGff($parm{'in_gff'}); 
	my %gidx = %{ &_index_glist($chr2gen) }; 
	
	my ($alnInfo) = &_readInAln( $parm{'in_aln'}, $opts{'useYN'} ); 
	my ($bn6Info) = &_readInBn6( $parm{'in_blast'} ); 
	$bn6Info = &_filter_eval( $bn6Info, $parm{'max_eval'} ) if ($parm{'max_eval'} >= 0); 
	$bn6Info = &_filter_cscore( $bn6Info, $parm{'min_cscore'}, $gen2loc ) if ($parm{'min_cscore'} > 0) ; 

	my ($prox_gene, $tand_gene); 


	if ( $parm{'only_direct'} ) {
		my %id2id; 
		for my $ar1 (@$bn6Info) {
			my $g1 = $ar1->{'k2v'}{'qseqid'}; 
			my $g2 = $ar1->{'k2v'}{'sseqid'}; 
			$id2id{$g1}{$g2} ++; 
			$id2id{$g2}{$g1} ++; 
		}

		for my $chrID (sort keys %$chr2gen) {
			my @gg = sort { $gidx{$a} <=> $gidx{$b} } map { $_->[0] } @{$chr2gen->{$chrID}}; 
			for (my $i=0; $i<@gg; $i++) {
				defined $id2id{$gg[$i]} or next; 
				my @tp = ($i); # proximal
				my @tt = ($i); # tandem
				for (my $j=$i+1; $j<@gg; $j++) {
					my $is_good = 0; 
					if ( $j-$i <= $parm{'max_proxN'} ) {
						defined $id2id{$gg[$i]}{$gg[$j]} and push(@tp, $j); 
						$is_good = 1; 
					}
					if ( $j-$i <= $parm{'max_tandN'} ) {
						defined $id2id{$gg[$i]}{$gg[$j]} and push(@tt, $j); 
						$is_good = 1; 
					}
					$is_good == 1 or last; 
				}
				scalar(@tp) > 1 and push(@{$prox_gene}, [@gg[@tp]]); 
				scalar(@tt) > 1 and push(@{$tand_gene}, [@gg[@tt]]); 
			}
		}
	} else {
		my ($id_grp) = &_group_id( $bn6Info ); 
		# [[grp0_gen1, grp0_gen2, ...], [grp1_gen1, grp2_gen2, ...], ...]
		$prox_gene = &_tandemGrp(
		 'id_grp'=>$id_grp, 
		 'max_distN'=>$parm{'max_proxN'}, 
		 'gen2loc'=>$gen2loc, 
		 'gindex'=>\%gidx, 
		); 
		$tand_gene = &_tandemGrp(
		 'id_grp'=>$id_grp, 
		 'max_distN'=>$parm{'max_tandN'}, 
		 'gen2loc'=>$gen2loc, 
		 'gindex'=>\%gidx, 
		); 
	}


	if ( defined $parm{'tand_list'} ) {
		my $otFh = &openFH($parm{'tand_list'}, '>'); 
		for my $ar1 (sort { $gidx{$a->[0]} <=> $gidx{$b->[0]} } @$tand_gene) {
			print {$otFh} join("\t", @$ar1)."\n"; 
		}
		close($otFh); 
	}
	if ( defined $parm{'prox_list'} ) {
		my $otFh = &openFH($parm{'prox_list'}, '>'); 
		for my $ar1 (sort { $gidx{$a->[0]} <=> $gidx{$b->[0]} } @$prox_gene) {
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


sub cvt_ctg2scf_byAGP {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	
	my %ctg2scf = %{ &fileSunhh::load_agpFile( $parm{'scf_agp'} ) }; 
	# %ctg2scf = ( $ctgID => [ [ctgS, ctgE, scfID, scfS, scfE, scfStr(+/-/?)], [], ... ] ) # This is sorted. 
	my ($alnInfo) = &_readInAln( $parm{'in_aln'} , $opts{'useYN'} ); 

	my %unknown_ctg; 

	defined $alnInfo->[0]{'text'} and print {$outFh} $alnInfo->[0]{'text'}; 
	for ( my $i=1; $i<@{$alnInfo}; $i++ ) {
		my ( $alnID, $n_pairs, $ctg1, $ctg2, $strand ) = @{$alnInfo->[$i]{'info'}}[0, 3, 4,5,6]; 
		$strand =~ m/^plus$/i and $strand = '+'; 
		$strand =~ m/^minus$/i and $strand = '-'; 
		$strand =~ m/^(\+|\-)$/i or &stopErr("[Err] Unknown strand [$strand]\n"); 
		my ($scf1_aR) = $ms_obj->switch_position( 'qry2ref' => \%ctg2scf, 'qryID' => $ctg1, 'qryPos' => 1, 'strand' => '+' ); 
		my ($scf2_aR) = $ms_obj->switch_position( 'qry2ref' => \%ctg2scf, 'qryID' => $ctg2, 'qryPos' => 1, 'strand' => $strand ); 
		unless ( defined $scf1_aR->[0] ) {
			defined $unknown_ctg{$ctg1} or do { &tsmsg("[Wrn] No contig information found for [$ctg1] in AGP file.\n"); $unknown_ctg{$ctg1} = 1; }; 
			$scf1_aR = [ $ctg1, 1, '+' ]; 
		}
		unless ( defined $scf2_aR->[0] ) {
			defined $unknown_ctg{$ctg2} or do { &tsmsg("[Wrn] No contig information found for [$ctg2] in AGP file.\n"); $unknown_ctg{$ctg2} = 1; }; 
			$scf2_aR = [ $ctg2, 1, '+' ]; 
		}

		$scf1_aR->[2] eq $scf2_aR->[2] or $strand =~ tr/+-/-+/; 
		my $rev = ( $scf1_aR->[2] eq '-' ) ? 1 : 0 ; 
		&out_blk( $alnInfo->[$i]{'text'}, $scf1_aR->[0], $scf2_aR->[0], $strand, $rev ); 

	}

}# cvt_ctg2scf_byAGP() 

sub out_blk {
	my ($blk_text, $new_ctg1, $new_ctg2, $strand, $reverse) = @_; 
	$strand eq '+' and $strand = 'plus'; 
	$strand eq '-' and $strand = 'minus'; 
	my @lines = split(/\n/, $blk_text); 
	my $n ; 

	if ($lines[0] =~ m!^(## Alignment \d+: score=\S+ e_value=\S+ )N=(\d+) ([^\&\s]+)\&([^\&\s]+) (plus|minus|X+)\s*$!) {
		my ( $part1, $n_pairs, $chr1, $chr2, $str ) 
		=  ( $1,     $2,       $3,    $4,    $5 ); 
		$n = $n_pairs; 
		my $new_line = "${part1}N=${n_pairs} $new_ctg1\&${new_ctg2} $strand"; 
		print {$outFh} "$new_line\n"; 
	} else {
		&stopErr("[Err] Unknown block title: $lines[0]\n"); 
	}

	unless ( $reverse ) {
		print {$outFh} join("\n", @lines[1 .. $#lines])."\n"; 
		return; 
	}

	for my $tline ( reverse(@lines[1 .. $#lines]) ) {
		$tline =~ m!^(\s*\d+)\-(\s*\d+):(\t\S+\t\S+\s*\S+.*)$! or &stopErr("[Err] Unknown block body: $tline\n"); 
		my ($i1, $i2, $part3) = ($1, $2, $3); 
		my $len2 = length($i2); 
		$i2 =~ s/\s//g; 
		$i2 = sprintf("%${len2}d", $n-$i2-1); 
		my $new_line = "${i1}-${i2}:$part3"; 
		print {$outFh} "$new_line\n"; 
	}
	
	return; 
}# out_blk() 

sub trim_block {
	my %parm = $ms_obj->_setHashFromArr(@_); 

	my $fh1 = &openFH( $parm{'trim_blkSE'}, '<' ); 
	my %toAdd; 
	my %toDelete; 
	my %toTrim; 
	my %processed; 
	while ( &wantLineC($fh1) ) {
		my @ta = &splitL("\t", $_); 
		if ( $ta[2] eq '.Remove' ) {
			$toDelete{$ta[0]}{$ta[1]} = 1; 
			$processed{'toDel'}{"$ta[0]\t$ta[1]"} = 0; 
		} elsif ( $ta[2] eq '.Add' ) {
			my @tb = split(/ : /, $ta[1]); 
			@tb == 2 or &stopErr("[Err] Bad input of field for .Add [$ta[1]]\n"); 
			my @tc = split(/ ; /, $tb[1]); 
			$tc[-1] =~ s/\s+$//; 
			my $txt = join("\t", @tc); 
			defined $processed{'toAdd'}{"$ta[0]\t$tb[0]\t$txt"} and do { &tsmsg("[Wrn] Skip repeated line : $_\n"); next; }; 
			push( @{$toAdd{$ta[0]}{$tb[0]}}, $txt ); 
			$processed{'toAdd'}{"$ta[0]\t$tb[0]\t$txt"} = 0; 
		} else {
			$toTrim{$ta[0]} = [ @ta[1,2] ]; 
			$processed{'toTrim'}{"$ta[0]\t$ta[1]\t$ta[2]"} = 0; 
		}
	}
	close($fh1); 

	my ($alnInfo) = &_readInAln( $parm{'in_aln'}, $opts{'useYN'} ); 

	defined $alnInfo->[0]{'text'} and print {$outFh} $alnInfo->[0]{'text'}; 

	for ( my $i=1; $i<@{$alnInfo}; $i++ ) {
		if (defined $toAdd{ $alnInfo->[$i]{'info'}[0] }) {
			my @ll = split(/\n/, $alnInfo->[$i]{'text'}); 
			my @new_ll = ($ll[0]); 
			my $new_cnt = -1; 

			for (my $j=1; $j<@ll; $j++) {
				$new_cnt ++; 
				my $gid1 = $alnInfo->[$i]{'pair'}[$j-1][0]; 
				$ll[$j] =~ m/\b${gid1}\b/ or &stopErr("[Err] geneID1 [$gid1] is not found in line : $ll[$j]\n"); 
				$ll[$j] =~ m/^(\s*\d+)\-(\s*\d+):(\t\S+\t\S+\s*\S+.*)$/ or &stopErr("[Err] Unknown block body: $ll[$j]\n"); 
				my ($i1, $i2, $part3) = ($1, $2, $3); 
				my $len2 = length($i2); 
				if ( $j == 1 and defined $toAdd{ $alnInfo->[$i]{'info'}[0] }{'.Top'} ) {
					for my $tl ( @{ $toAdd{ $alnInfo->[$i]{'info'}[0] }{'.Top'} } ) {
						$i2 = sprintf("%${len2}d", $new_cnt); 
						push(@new_ll, "${i1}-${i2}:\t$tl"); 
						$new_cnt++; 
						$processed{'toAdd'}{"$alnInfo->[$i]{'info'}[0]\t.Top\t$tl"} = 1; 
					}
				}
				$i2 = sprintf("%${len2}d", $new_cnt); 
				push(@new_ll, "${i1}-${i2}:$part3"); 
				if ( defined $toAdd{ $alnInfo->[$i]{'info'}[0] }{$gid1} ) {
					for my $tl ( @{ $toAdd{ $alnInfo->[$i]{'info'}[0] }{$gid1} } ) {
						$new_cnt ++; 
						$i2 = sprintf("%${len2}d", $new_cnt); 
						push(@new_ll, "${i1}-${i2}:\t$tl"); 
						$processed{'toAdd'}{"$alnInfo->[$i]{'info'}[0]\t$gid1\t$tl"} = 1; 
					}
				}
			}
			$alnInfo->[$i]{'text'} = join("\n", @new_ll)."\n"; 
			&_update_alnInfo_byText( $alnInfo->[$i] ); 
		}
		if (defined $toDelete{ $alnInfo->[$i]{'info'}[0] }) {
			my @ll = split(/\n/, $alnInfo->[$i]{'text'}); 
			my @new_ll = ($ll[0]); 
			my $new_cnt = -1; 
			for (my $j=1; $j<@ll; $j++) {
				my $gid1 = $alnInfo->[$i]{'pair'}[$j-1][0]; 
				if ( defined $toDelete{ $alnInfo->[$i]{'info'}[0] }{$gid1} ) {
					$processed{'toDel'}{"$alnInfo->[$i]{'info'}[0]\t$gid1"} = 1; 
					next; 
				}
				$new_cnt ++; 
				$ll[$j] =~ m/\b${gid1}\b/ or &stopErr("[Err] geneID1 [$gid1] is not found in line : $ll[$j]\n"); 
				$ll[$j] =~ m/^(\s*\d+)\-(\s*\d+):(\t\S+\t\S+\s*\S+.*)$/ or &stopErr("[Err] Unknown block body: $ll[$j]\n"); 
				my ($i1, $i2, $part3) = ($1, $2, $3); 
				my $len2 = length($i2); 
				$i2 = sprintf("%${len2}d", $new_cnt); 
				push(@new_ll, "${i1}-${i2}:$part3"); 
			}
			$alnInfo->[$i]{'text'} = join("\n", @new_ll)."\n"; 
			&_update_alnInfo_byText( $alnInfo->[$i] ); 
		}


		defined $toTrim{ $alnInfo->[$i]{'info'}[0] } or do { print {$outFh} $alnInfo->[$i]{'text'} ; next; }; 
		my $blkID = $alnInfo->[$i]{'info'}[0]; 
		my ( $sGID, $eGID ) = @{$toTrim{ $alnInfo->[$i]{'info'}[0] }}; 
		$processed{'toTrim'}{"$alnInfo->[$i]{'info'}[0]\t$sGID\t$eGID"} = 1; 
		my ($sj, $ej); 
		my @ll = split(/\n/, $alnInfo->[$i]{'text'}); 
		for (my $j=1; $j<@ll; $j++) {
			my $gid1 = $alnInfo->[$i]{'pair'}[$j-1][0]; 
			my $gid2 = $alnInfo->[$i]{'pair'}[$j-1][1]; 
			$ll[$j] =~ m/\b${gid1}\b/ or &stopErr("[Err] geneID1 [$gid1] is not found in line : $ll[$j]\n"); 
			$gid1 eq $sGID and $sj = $j; 
			$gid1 eq $eGID and $ej = $j; 
		}
		defined $sj or do { &tsmsg("[Wrn] geneID1 [$sGID] is not found in block [$blkID]\n"); $sj = 1; }; 
		defined $ej or do { &tsmsg("[Wrn] geneID1 [$eGID] is not found in block [$blkID]\n"); $ej = $#ll; }; 
		$sj > $ej and ($sj, $ej) = ($ej, $sj); 
		my $newN = $ej-$sj+1; 
		$ll[0] =~ m!^(## Alignment \d+: score=\S+ e_value=\S+ )N=(\d+) ([^\&\s]+\&[^\&\s]+ (?i:plus|minus|X+))\s*$! or &stopErr("[Err] Unknown block line : $ll[0]\n"); 
		print {$outFh} "$1N=${newN} $3\n"; 
		for my $tline ( @ll[$sj .. $ej] ) {
			$tline =~ m!^(\s*\d+)\-(\s*\d+):(\t\S+\t\S+\s*\S+.*)$! or &stopErr("[Err] Unknown block body: $tline\n"); 
			my ($i1, $i2, $part3) = ($1, $2, $3);
			my $len2 = length($i2);
			$i2 =~ s/\s//g;
			$i2 = sprintf("%${len2}d", $i2-$sj+1); 
			my $new_line = "${i1}-${i2}:$part3"; 
			print {$outFh} "$new_line\n";
		}
	}

	# Check if all processed. 
	for my $k1 (sort keys %processed) {
		for my $k2 (sort keys %{$processed{$k1}}) {
			$processed{$k1}{$k2} == 0 and &tsmsg("[Wrn] Infor for [$k1][$k2] not processed.\n"); 
		}
	}

	return(); 
}# trim_block() 

sub _update_alnInfo_byText {
	my ($hr) = @_; 
	my @ll = split(/\n/, $hr->{'text'}); 
	$hr->{'info'}[3] == $#ll or do { &tsmsg("[Msg] Update Num_pairs to $#ll in Block [@{$hr->{'info'}}]\n"); $hr->{'info'}[3] = $#ll; }; 
	$ll[0] =~ m!^(## Alignment \d+: score=\S+ e_value=\S+ )N=(\d+) ([^\&\s]+\&[^\&\s]+ (?i:plus|minus|X+))\s*$! or &stopErr("[Err] Unknown block line : $ll[0]\n");
	$ll[0] = "$1N=$#ll $3"; 
	$hr->{'text'} = join("\n", @ll)."\n"; 
	my @new_pair; 
	for my $tline (@ll) {
		if ( $tline =~  m!^\s*(\d+)\-\s*(\d+):\t(\S+)\t(\S+)\s*(\S+)(?:\t(\S+)\t(\S+)(?:\t(\S+))?)?$! ) {
			# This is the normal format or the KaKs output of MCscanX. 
			my ($alnID, $alnID_id, $gid1, $gid2, $eval, $tka, $tks, $tw) 
			= 
			(   $1,     $2,        $3,    $4,    $5,    $6,   $7,   $8);
			$tka //= ''; $tks //= ''; $tw //= '';
			push(@new_pair, [$gid1, $gid2, $eval, $tka, $tks, $tw]);
		} elsif ( $tline =~ m!^\s*(\d+)\-\s*(\d+):\t(\S+)\t(\S+)\s*(\S+)(?:\t(\S+)\t(\S+)(?:\t(\S+)\t(\S+))?)?$! ) {
			# This is the format of KsKa addition from my follow_msc perl script. 
			# #  0-  0:  Cma_000007      Cma_000973        2e-57 0.5258  0.0154  0.5387  0.0148
			my ($alnID, $alnID_id, $gid1, $gid2, $eval, $tks, $tka, $t_ngKs, $t_ngKa)
			= 
			(   $1,     $2,        $3,    $4,    $5,    $6,   $7,   $8,      $9);
			my $tw;
			if ( defined $t_ngKs and !$opts{'useYN'} ) {
				$tks = $t_ngKs; 
				$tka = $t_ngKa; 
			}
			$tks < 0 and $tks = 'nan'; 
			if ( defined $tks and $tks ne 'nan') {
				$tw = ( $tks > 0 ) ? $tka/$tks : 'nan';
			}
			$tks eq 'nan' and $tw = 'nan';
			$tka //= ''; $tks //= ''; $tw //= '';
			push(@new_pair, [$gid1, $gid2, $eval, $tka, $tks, $tw]);
		}
	}
	$hr->{'pair'} = [@new_pair]; 
	return; 
}# _update_alnInfo_byText() 

sub mcs_blkByList {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	
	my ($alnInfo) = &_readInAln($parm{'in_aln'}, $opts{'useYN'}); 
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

sub filter_block {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	$parm{'min_blkPair'} //= 2; 
	
	my ($alnInfo) = &_readInAln($parm{'in_aln'}, $opts{'useYN'}); 

	my @good_idx; 
	for (my $i=0; $i<@$alnInfo; $i++) {
		if ($i == 0 and !(defined $alnInfo->[$i]{'info'})) {
			push(@good_idx, $i); 
			print {$outFh} $alnInfo->[$i]{'text'}; 
			next; 
		}
		$alnInfo->[$i]{'info'}[3] >= $parm{'min_blkPair'} or next; 
		push(@good_idx, $i); 
		print {$outFh} $alnInfo->[$i]{'text'}; 
	}

	return; 
}# filter_block ()

sub mcs_aln2pair {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	$parm{'useYN'} //= 0; 

	my ($alnInfo) = &_readInAln($parm{'in_aln'}, $parm{'useYN'}); 
	defined $alnInfo->[0]{'info'} or shift(@{$alnInfo}); 

	print join("\t", qw/BlkID orderID Gene1 Gene2 Ka Ks w/)."\n"; 
	for (my $i=0; $i<@{$alnInfo}; $i++) {
		my (@gen1, @gen2, @ka, @ks, @w); 
		for (my $j=0; $j<@{$alnInfo->[$i]{'pair'}}; $j++) {
			my $ar1 = $alnInfo->[$i]{'pair'}[$j]; 
			print join("\t", $alnInfo->[$i]{'info'}[0], $j, $ar1->[0], $ar1->[1], $ar1->[3], $ar1->[4], $ar1->[5])."\n"; 
		}
	}

	return; 
} # mcs_aln2pair

sub mcs_table2aln {
  my %parm = $ms_obj->_setHashFromArr(@_);

  open F,'<',"$parm{'in_table'}" or die;
  while (<F>) {
    chomp;
    my @ta=split(/\t/, $_);
    if ( $ta[0] eq 'BlkID' ) {
      next;
    }
    my ($alnID, $alnScore, $alnEval, $alnNum, $chr1, $chr2, $str) = @ta[0, 8, 9, 10, 1, 4, 7];
    $str = ($str eq "-") ? "minus" : "plus" ;
    print STDOUT "## Alignment $alnID: score=$alnScore e_value=$alnEval N=$alnNum $chr1\&$chr2 $str\n";
    my @gene1 = split(/,/, $ta[11]);
    my @gene2 = split(/,/, $ta[12]);
    for (my $i=0; $i<@gene1; $i++) {
      print STDOUT join("\t", "$alnID\-$i:", $gene1[$i], $gene2[$i], "0.0")."\n";
    }
  }
  close F;
}# mcs_table2aln()

sub msc_aln2table {
	my %parm = $ms_obj->_setHashFromArr(@_); 
	
	my ($chr2gen, $gen2loc)  = &_readInGff($parm{'in_gff'}); 
	my ($alnInfo) = &_readInAln($parm{'in_aln'}, $opts{'useYN'}); 
	defined $alnInfo->[0]{'info'} or shift(@{$alnInfo}); 
	
	print {$outFh} join("\t", qw/BlkID Chrom1 Start1 End1 Chrom2 Start2 End2 Strand AlnScore AlnEvalue AlnNumber Gene1 Gene2 Ka Ks KaKs AVG_Ks Med_Ks UsedKs AVG_Ka Med_Ka UsedKa AVG_KaKs Med_KaKs UsedKaKs INSavg_Ks INSmed_Ks INSused_Ks/)."\n"; 
	for (my $i=0; $i<@{$alnInfo}; $i++) {
		my (@gen1, @gen2, @ka, @ks, @w); 
		for my $ar1 (@{$alnInfo->[$i]{'pair'}}) {
			push(@gen1, $ar1->[0]); 
			push(@gen2, $ar1->[1]); 
			push(@ka, $ar1->[3]); 
			push(@ks, $ar1->[4]); 
			push(@w,  $ar1->[5]); 
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
		 &_avg(@ks), &_avg(@ka), &_avg(@w), &_avgINS(@ks) 
		)."\n"; 
	}
}# msc_aln2table() 

sub msc_aln2list{
	my %parm = $ms_obj->_setHashFromArr(@_); 
	$parm{'addChr'} //= 0; 
	$parm{'srt_by'} //= 'min'; 
	$parm{'raw_order'} //= 0; 
	
	my ($chr2gen, $gen2loc) = &_readInGff($parm{'in_gff'}); 
	my ($chr2gen_tgt, $gen2loc_tgt) = ($chr2gen, $gen2loc); 
	if (defined $parm{'tgt_gff'} and $parm{'tgt_gff'} ne '') {
		($chr2gen_tgt, $gen2loc_tgt) = &_readInGff($parm{'tgt_gff'}); 
	}

	my ($alnInfo) = &_readInAln($parm{'in_aln'}, $opts{'useYN'}); 
	defined $alnInfo->[0]{'info'} or shift(@{$alnInfo}); 

	print {$outFh} join("\t", qw/Chromosome IntraDep InterDep InterTax Pivot PivotStart Blocks/)."\n"; 
	for my $chrID ( sort keys %$chr2gen ) {
		# ([genID, chrS], [], ...)
		my @glist; 
		unless ( $parm{'raw_order'} ) {
			@glist = map { [ $_->[0], $_->[1] ] } sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @{$chr2gen->{$chrID}}; 
		} else {
			@glist = map { [ $_->[0], $_->[1] ] } @{$chr2gen->{$chrID}}; 
		}
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
				$use_tax{$tax[$nb]} //= 1; 

defined $parm{'gindex'}{ $gpair[$k][$na] } or &stopErr("[Err] No gene location found for [$gpair[$k][$na]]\n"); 
my $needIdx = $parm{'gindex'}{ $gpair[$k][$na] }; 
# &tsmsg("[Msg] j=$j needIdx=$needIdx $gpair[$k][$na] $gpair[$k][$nb]\n"); 
if ($j > $needIdx) {
	&tsmsg("[Err] @{$tr->{'info'}}\n"); 
	&tsmsg("[Err] minIdx=$minIdx maxIdx=$maxIdx\n"); 
	&tsmsg("[Err] gene [k=$k na=$na $gpair[$k][$na]][nb=$nb $gpair[$k][$nb]]\n"); 
	&tsmsg("[Err] depth=$depth[$j]{$tax[$nb]}; use_tax=$use_tax{$tax[$nb]}; tax_nb=$tax[$nb]\n"); 
	&tsmsg("[Err] The location of gene [$gpair[$k][$na]] might be different from .aln file. [j=$j needIdx=$needIdx]\n"); 
	&stopErr("[Err] You may try to add -raw_order parameter.\n"); 
}

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

sub _readInAln {
	return( &mcsSunhh::_readInAln( @_ ) ); 
}# _readInAln() 

# Return : (\@alnInfo)
#  @alnInfo : 
#    [idx_num]{'info'} => [alnID, score, evalue, Num_pairs, chrID_1, chrID_2, strand(plus|minus)]
#    [idx_num]{'pair'} => [ [genID_1, genID_2, pair_evalue, Ka, Ks, Ka/Ks], [], ... ]
#    [idx_num]{'text'} => $text_of_current_block
sub _readInAln_1{
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
			$aln_idx > 0 or &stopErr("Too early to line: $_\n"); 
			my ($alnID, $alnID_id, $gid1, $gid2, $eval, $tka, $tks, $tw) 
			= 
			   ($1,     $2,        $3,    $4,    $5,    $6,   $7,   $8); 
			$tka //= ''; $tks //= ''; $tw //= ''; 
			$alnID == $alnInfo[$aln_idx]{'info'}[0] or &stopErr("[Err] line_alnID=$alnID not fitting upper level (alnID=$alnInfo[$aln_idx]{'info'}[0]).\n"); 
			push(@{$alnInfo[$aln_idx]{'pair'}}, [$gid1, $gid2, $eval, $tka, $tks, $tw]); 
			$alnInfo[$aln_idx]{'text'} .= "$_\n"; 
		} elsif ( m!^\s*(\d+)\-\s*(\d+):\t(\S+)\t(\S+)\s*(\S+)(?:\t(\S+)\t(\S+)(?:\t(\S+)\t(\S+))?)?$! ) { 
			# #  0-  0:  Cma_000007      Cma_000973        2e-57 0.5258  0.0154  0.5387  0.0148
			#                                                    yn_Ks   yn_Ka   ng_Ks   ng_Ka
			# This is to fit result from haibao tang's python ks calculation. 
			$aln_idx > 0 or &stopErr("Too early to line: $_\n"); 
			my ($alnID, $alnID_id, $gid1, $gid2, $eval, $tks, $tka, $t_ngKs, $t_ngKa) 
			= 
			(   $1,     $2,        $3,    $4,    $5,    $6,   $7,   $8,      $9); 
			my $tw; 
			if ( !$opts{'useYN'} and defined $t_ngKs ) {
				$tks = $t_ngKs; 
				$tka = $t_ngKa; 
			}
			$tks < 0 and $tks = 'nan'; 
			if ( defined $tks and $tks ne 'nan') {
				$tw = ( $tks > 0 ) ? $tka/$tks : 'nan'; 
			}
			$tks eq 'nan' and $tw = 'nan'; 
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
}# _readInAln_1

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

sub _avgINS {
	my ($avg, $med) = ('nan', 'nan'); 
	my @unn = grep { defined $_ and $_ !~ m!^(|\-?nan|\-?inf)$!i and $_ >= 0 } @_; 
	my $validNum = scalar(@unn); 
	if ( $validNum > 0 ) {
		my $ic = $ms_obj->ins_calc(\@unn, 0); 
		$avg = $ic->{'interval_mean'}; 
		$med = $ic->{'interval_median'}; 
		$validNum = $ic->{'interval_cnt'}; 
	}
	return ($avg, $med, $validNum); 
}# _avgINS() 

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
		my @tb; 
		if ( $opts{'raw_order'} ) {
			@tb = @{$chr2gen->{$chrID}}; 
		} else {
			@tb = sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @{$chr2gen->{$chrID}}; 
		}
		for my $ar1 ( @tb ) {
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
		             || ( ( $opts{'raw_order'} ) ? -1 
                                             : ($parm{'gen2loc'}{$a}[1] <=> $parm{'gen2loc'}{$b}[1]
                                               || $parm{'gen2loc'}{$a}[2] <=> $parm{'gen2loc'}{$b}[2] 
                                               ) 
                    )
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

