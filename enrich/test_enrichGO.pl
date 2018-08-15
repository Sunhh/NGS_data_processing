#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"inSubList:s", 
	  "colN:i", # 0 
	  "skipHlineN:i", # 1 
	"padj:f", 
	"opref:s", 

	"bgGOtab:s", 
	"oboTab:s", 
); 
$opts{'colN'}       //= 0; 
$opts{'skipHlineN'} //= 0; 
$opts{'padj'}       //= 0.05; 
$opts{'opref'}      //= 'out'; 

my $help_txt = <<"H1"; 
################################################################################
# perl $0 background_pwy subset_geneList [padj_cutoff_0.05]
# perl $0  -inSubList  subset_geneList  -bgGOtab wm97pbV2ID_feiID_clean.GO_bg.tab  -oboTab gene_ontology_edit.obo.2018-05-01.tab 
#
# -help 
# 
# -colN        [$opts{'colN'}] Column number of gene ID. 
# -skipHlineN  [$opts{'skipHlineN'}] Number of lines to be ignored in subset list. 
# 
# -padj        [$opts{'padj'}] FDR threshold for output. Using '<' instead of '<='
#
# -inSubList   [filename] Required. A list of subset genes. 
#
# -opref       ['$opts{'opref'}'] Prefix of output files. .tab; 
# 
# -bgGOtab     [filename] Required. Come from commands : 
#                perl split_jnGO.pl -skipHlineN 1 -in wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc > wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc.2col
#                cat wm97pbV2ID_evm_feiID_clean.byChr.gff3.jnLoc | tail -n +2 | deal_table.pl -column 0 | deal_table.pl -kSrch_idx wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc.2col -kSrch_drop | deal_table.pl -label_mark NA | deal_table.pl -column 1,0,0 > wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc.2col.addFakeRest
#                cat wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc.2col wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc.2col.addFakeRest > wm97pbV2ID_feiID_clean.GO_bg.tab
#                # Example : 
#                eleID                   GO_ID           GO_descTxt
#                Cla97C01G000010.1       GO:0016021      integral component of membrane
#                Cla97C01G000030.1       GO:0004222      metalloendopeptidase activity
#                Cla97C01G000030.1       GO:0006508      proteolysis
#                Cla97C01G0000XX.1       NA              NA
#
#
# -oboTab      [filename] Required. Come from command : 
#                gzip -cd gene_ontology_edit.obo.2018-05-01.gz | perl cnvt_GOobo_to_tab.pl > gene_ontology_edit.obo.2018-05-01.tab
################################################################################
H1

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
for my $k1 (qw/inSubList bgGOtab oboTab/) {
	defined $opts{$k1} or do { &tsmsg("[Err] Lack -$k1\n"); &LogInforSunhh::usage($help_txt); }; 
}

my (%oboInfo, %goBg, %subInfo); 
&load_oboTab(); 
&load_goBg(); 
&load_subset(); 
&run_fisherET(); 

sub run_fisherET {
	# Run fisher's exact test. 
	my $wrkDir = &fileSunhh::new_tmp_dir('create'=>1); 
	my $fh_oTab = &openFH("$wrkDir/input", '>'); 
	print {$fh_oTab} join("\t", qw/GO_ID n11 n12 n21 n22/)."\n"; 
	for my $goid (keys %{$subInfo{'goid_cntInSubset'}}) {
		# In R - fisher.test() : 
		### total ele in genome : $goBg{'eleNum'}
		### GO    ele in genome : $goBg{'goid_cntInGenome'}{$goid}
		### total ele in subset : $subInfo{'ele_cntInSubset'}
		### GO    ele in subset : $subInfo{'goid_cntInSubset'}{$goid}
		my $n11 = $subInfo{'goid_cntInSubset'}{$goid} ; # in GO & in subset
		my $n12 = $subInfo{'ele_cntInSubset'} - $subInfo{'goid_cntInSubset'}{$goid}; # not GO & in subset
		my $n21 = $goBg{'goid_cntInGenome'}{$goid} - $subInfo{'goid_cntInSubset'}{$goid}; # in GO & not subset
		my $n22 = $goBg{'eleNum'} - $n11 - $n12 - $n21; # not GO & not subset
		print {$fh_oTab} join("\t", $goid, $n11, $n12, $n21, $n22)."\n"; 
	}
	close($fh_oTab); 
	my $fh_r = &openFH("$wrkDir/c.r", '>'); 
	print {$fh_r} <<"R1"; 
library(stats)
tbl     <- read.table( '$wrkDir/input', header=TRUE, sep="\\t", stringsAsFactors=FALSE, row.names=1 )
rawP    <- apply( tbl, MARGIN=1, FUN= function(x) { mat<-matrix( x[1:4], nrow=2, byrow=TRUE ); fisher.test(mat, alternative="greater")\$p.value } )
adjP    <- p.adjust( rawP, method= 'BH' ) # BH == FDR; 
ooo     <- data.frame( goid=rownames(tbl), rawP=rawP, FDR=adjP )
ooo     <- ooo[ sort(ooo\$FDR, decreasing=T, index.return=TRUE)\$ix, ]
write.table( ooo, file='$wrkDir/output', append=FALSE, quote=FALSE, sep="\\t", row.names=FALSE, col.names=FALSE )

R1
	&exeCmd_1cmd("Rscript $wrkDir/c.r"); 

	my $fh_i = &openFH("$wrkDir/output", '<'); 
	my @out_lines; 
	while (<$fh_i>) {
		chomp; 
		m!^\s*$! and next; 
		my @ta = &splitL("\t", $_); 
		my ($goid, $rawP, $fdr)  = @ta; 
		$fdr < $opts{'padj'} or next; 
		push(@out_lines, 
		[
			$goid, $fdr, $rawP, 
			$oboInfo{'goid2type'}{$goid}, 
			$oboInfo{'goid2name'}{$goid}, 
			$oboInfo{'goid2def'}{$goid}, 
			join(";", @{$subInfo{'goid2ele'}{$goid}}), 
			join("_", 
				$subInfo{'goid_cntInSubset'}{$goid}, 
				$subInfo{'ele_cntInSubset'}, 
				$goBg{'goid_cntInGenome'}{$goid}, 
				$goBg{'eleNum'}
			)
		]); 
	}
	close($fh_i); 
	my $fh_o = &openFH("$opts{'opref'}.tab", '>'); 
	@out_lines > 0 and print {$fh_o} join("\t", qw/GO_ID FDR rawP GO_Type GO_name GO_def GO_eleID_sub sub_goN_totalN_bg_goN_totalN/)."\n"; 
	for (sort { $a->[3] cmp $b->[3] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @out_lines) {
		print {$fh_o} join("\t", @$_)."\n"; 
	}
	close($fh_o); 
	&fileSunhh::_rmtree($wrkDir); 
	return; 
}# run_fisherET() 


####################################################################################################
#  sub routines. 
####################################################################################################
sub load_oboTab {
	# %oboInfo: goid2name, goid2def, goid2type, altid2id, gotype2id; 
	my $fh = &openFH($opts{'oboTab'}, '<'); 
	while (<$fh>) {
		chomp; 
		my @ta = &splitL("\t", $_); 
		my ($goid, $altid, $type, $goname, $godef) = @ta; 
		$goid eq 'GO_ID' and next; 
		$type eq 'ET' and next; 
		$oboInfo{'goid2name'}{$goid} = $goname; 
		$oboInfo{'goid2def'}{$goid}  = $godef; 
		$oboInfo{'goid2type'}{$goid} = $type; 
		$oboInfo{'altid2id'}{$goid}  = $goid; 
		push(@{$oboInfo{'gotype2id'}{$type}}, $goid); 
		unless ($ta[1] =~ m!^\s*NA\s*$!i) {
			for my $t1 (split(/;/, $altid)) {
				$oboInfo{'altid2id'}{$t1} = $goid; 
			}
		}
	}
	close($fh); 
	return; 
}# load_oboTab
sub load_goBg {
	# %goBg: eleNum , ele2goid, goid2ele
	my $fh = &openFH($opts{'bgGOtab'}, '<'); 
	my %havEleID; 
	my %pair_ele_go; 
	while (<$fh>) {
		chomp; 
		my @ta = &splitL("\t", $_); 
		my ($eleID, $goid, $go_desc) = @ta; 
		$eleID =~ m!^(eleID|geneID|mrnaID)$!i and next; 
		$havEleID{$eleID} = 1; 
		$goid =~ m!^\s*NA\s*$!i and next; 
		if ( !( defined $oboInfo{'altid2id'}{$goid} ) ) {
			&tsmsg("[Wrn] No information for GO_ID [$goid]\n"); 
			for my $t0 (qw/goid2type goid2name goid2def/) {
				$oboInfo{$t0}{$goid} = 'NA'; 
			}
			$oboInfo{'altid2id'}{$goid} = $goid; 
		} else {
			$goid = $oboInfo{'altid2id'}{$goid}; 
		}
		my $tk = "$eleID\t$goid"; 
		defined $pair_ele_go{$tk} and next; 
		$pair_ele_go{$tk} = 1; 
		push(@{$goBg{'ele2goid'}{$eleID}}, $goid); 
	}
	close($fh); 
	$goBg{'eleNum'} = scalar(keys %havEleID); 
	for my $eleID (keys %{$goBg{'ele2goid'}}) {
		for my $goid (@{$goBg{'ele2goid'}{$eleID}}) {
			push(@{$goBg{'goid2ele'}{$goid}}, $eleID); 
		}
	}
	for my $goid (keys %{$goBg{'goid2ele'}}) {
		$goBg{'goid_cntInGenome'}{$goid} = scalar( @{$goBg{'goid2ele'}{$goid}} ); 
	}
	return; 
}# load_goBg
sub load_subset {
	# %subInfo: goid_cntInSubset , ele_cntInSubset ; 
	my $fh = &openFH($opts{'inSubList'}, '<'); 
	my %havEleID; 
	while (<$fh>) {
		$. <= $opts{'skipHlineN'} and next; 
		chomp; 
		my @ta = &splitL("\t", $_); 
		( defined $ta[$opts{'colN'}] and $ta[$opts{'colN'}] !~ m!^\s*(NA|)\s*$!i) or next; 
		my $eleID = $ta[$opts{'colN'}]; 
		defined $havEleID{$eleID} and next; 
		$havEleID{$eleID} = 1; 
		if (defined $goBg{'ele2goid'}{$eleID}) {
			for my $goid (@{$goBg{'ele2goid'}{$eleID}}) {
				$subInfo{'goid_cntInSubset'}{$goid} ++; 
				push(@{$subInfo{'ele2goid'}{$eleID}}, $goid); 
				push(@{$subInfo{'goid2ele'}{$goid}}, $eleID); 
			}
		}
	}
	close($fh); 
	$subInfo{'ele_cntInSubset'} = scalar(keys %havEleID); 
	scalar(keys %{$subInfo{'goid_cntInSubset'}}) == 0 and die "ffffff\n"; 
	return; 
}# load_subset


