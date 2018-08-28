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

	"enrichType:s", # Default 'GO'; 
	"bgGOtab:s", 
	"bgTab:s", 
	"oboTab:s", 
); 

my %gg;     # All global parameters. 
my %grpBg;  # group information in background; 
my %grpSub; # group information in subset; 

&set_Glob(); 
### For GO enrichment analysis ; 
my (%oboInfo); 
if ( $gg{'enrichType'} eq 'go' ) {
	&load_oboTab(); 
	&load_goBg(); 
	&load_subset(); 
	&run_fisherET(); 
} elsif ( $gg{'enrichType'} eq 'kegg' ) {
	&load_goBg(); 
	&load_subset(); 
	&run_fisherET(); 
} else {
	&tsmsg("[Err] Unsupported -enrichType $gg{'enrichType'}\n"); 
}


####################################################################################################
#  sub routines. 
####################################################################################################
sub set_Glob {
	$gg{'colN'}                  = 0; 
	$gg{'skipHlineN'}            = 0; 
	$gg{'padj'}                  = 0.05; 
	$gg{'opref'}                 = 'out'; 
	$gg{'enrichType'}            = 'GO'; 
	$gg{'htxt'}                  = <<"H1"; 
################################################################################
# perl $0  -inSubList  subset_geneList  -bgGOtab wm97pbV2ID_feiID_clean.GO_bg.tab  -oboTab gene_ontology_edit.obo.2018-05-01.tab 
#
# -help 
# 
# -colN        [$gg{'colN'}] Column number of gene ID. 
# -skipHlineN  [$gg{'skipHlineN'}] Number of lines to be ignored in subset list. 
# 
# -padj        [$gg{'padj'}] FDR threshold for output. Using '<' instead of '<='
#
# -inSubList   [filename] Required. A list of subset genes. 
#
# -opref       ['$gg{'opref'}'] Prefix of output files. .tab; 
#
# -enrichType  ['$gg{'enrichType'}'] Enrichment type. 
# 
#   enrichType = GO : 
#    -bgGOtab     [filename] Required. Come from commands : 
#                   perl split_jnGO.pl -skipHlineN 1 -in wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc > wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc.2col
#                   cat wm97pbV2ID_evm_feiID_clean.byChr.gff3.jnLoc | tail -n +2 | deal_table.pl -column 0 | deal_table.pl -kSrch_idx wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc.2col -kSrch_drop | deal_table.pl -label_mark NA | deal_table.pl -column 1,0,0 > wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc.2col.addFakeRest
#                   cat wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc.2col wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc.2col.addFakeRest > wm97pbV2ID_feiID_clean.GO_bg.tab
#                   # Example : 
#                   eleID                   GO_ID           GO_descTxt
#                   Cla97C01G000010.1       GO:0016021      integral component of membrane
#                   Cla97C01G000030.1       GO:0004222      metalloendopeptidase activity
#                   Cla97C01G000030.1       GO:0006508      proteolysis
#                   Cla97C01G0000XX.1       NA              NA
#
#
#    -oboTab      [filename] Required. Come from command : 
#                   gzip -cd gene_ontology_edit.obo.2018-05-01.gz | perl cnvt_GOobo_to_tab.pl > gene_ontology_edit.obo.2018-05-01.tab
#   enrichType = KEGG : 
#    -bgTab       [filename] Required. 
################################################################################
H1
	$opts{'help'} and &LogInforSunhh::usage($gg{'htxt'}); 
	for my $optK (qw/colN skipHlineN padj opref enrichType/) {
		defined $opts{$optK} and $gg{$optK} = $opts{$optK}; 
	}
	$gg{'enrichType'} = lc($gg{'enrichType'}); 
	for my $k1 (qw/inSubList/) {
		defined $opts{$k1} or do { &tsmsg("[Err] Lack -$k1\n"); &LogInforSunhh::usage($gg{'htxt'}); };
		$gg{$k1} = $opts{$k1}; 
	}
	if ($gg{'enrichType'} eq 'go') {
		for my $k2 (qw/bgGOtab oboTab/) {
			defined $opts{$k2} or do { &tsmsg("[Err] Lack -$k2\n"); &LogInforSunhh::usage($gg{'htxt'}); };
			$gg{$k2} = $opts{$k2}; 
		}
	} elsif ($gg{'enrichType'} eq 'kegg') {
		if ( defined $opts{'bgTab'} ) {
			$gg{'bgTab'} = $opts{'bgTab'}; 
		} elsif ( defined $opts{'bgGOtab'} ) {
			$gg{'bgTab'} = $opts{'bgGOtab'}; 
		} else {
			&tsmsg("[Err] Lack -bgTab\n"); 
			&LogInforSunhh::usage($gg{'htxt'}); 
		}
	} else {
		&stopErr("[Err] Unknown -enrichType [$gg{'enrichType'}]!\n"); 
	}
}# set_Glob() 

sub load_oboTab {
	# Fill %oboInfo: goid2name, goid2def, goid2type, altid2id, gotype2id; 
	my $fh = &openFH($gg{'oboTab'}, '<'); # Output of : gzip -cd gene_ontology_edit.obo.2018-05-01.gz | perl cnvt_GOobo_to_tab.pl > gene_ontology_edit.obo.2018-05-01.tab
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
	# I want to change goBg to a general structure of %grpBg
	#   $grpBg{'ele2grpID'} = { eleID=>[grpID1, grpID2, ...], ... } 
	#   $grpBg{'grpID2ele'} = { grpID=>[eleID1, eleID2, ...], ... }
	#   $grpBg{'eleID'}     = { eleID=>1 }
	#   $grpBg{'eleNum'}    = scalar(keys %{$grpBg{'ele2grpID'}}) 
	#   $grpBg{'grpSize_inGenome'} = { grpID1=>NumberOfEleInGrpInGenomeBg, grpID2=> , ... }
	#   $grpBg{'grpID2def'} = { grpID=>'group_description' }
	my $fh; 
	if ($gg{'enrichType'} eq 'go') {
		$fh = &openFH($gg{'bgGOtab'}, '<'); # Storing GO background information 
	} elsif ($gg{'enrichType'} eq 'kegg') {
		$fh = &openFH($gg{'bgTab'}, '<'); 
	} else {
		&tsmsg("[Err] Bad -enrichType [$gg{'enrichType'}]\n"); 
	}
	my %pair_ele_go; 
	while (<$fh>) {
		chomp; 
		my @ta = &splitL("\t", $_); 
		my ($eleID, $goid, $go_desc) = @ta; 
		$eleID =~ m!^(eleID|geneID|mrnaID)$!i and next; 
		$grpBg{'eleID'}{$eleID} = 1; 
		# Ignore GO ID with 'NA'; 
		$goid =~ s!^\s+|\s+$!!g; 
		if ( $goid =~ m!^NA$!i ) {
			next; 
		}
		if ( $gg{'enrichType'} eq 'go' ) {
			# Convert current goid into final GO ID. 
			if ( !( defined $oboInfo{'altid2id'}{$goid} ) ) {
				&tsmsg("[Wrn] No information for GO_ID [$goid]\n"); 
				for my $t0 (qw/goid2type goid2name goid2def/) {
					$oboInfo{$t0}{$goid} = 'NA'; 
				}
				$oboInfo{'altid2id'}{$goid} = $goid; 
			} else {
				$goid = $oboInfo{'altid2id'}{$goid}; 
			}
		} elsif ( $gg{'enrichType'} eq 'kegg' ) {
			; 
		} else {
			&tsmsg("[Err] Unsupported -enrichType '$gg{'enrichType'}'\n"); 
		}
		# Record group ID and its definition. 
		$grpBg{'grpID2def'}{$goid} //= $go_desc; 
		$grpBg{'grpID2def'}{$goid} eq $go_desc or &stopErr("[Err] Different definitions for grpID [$goid]: '$grpBg{'grpID2def'}{$goid}' VS. '$go_desc'\n"); 
		# Skip repeated GO ID; 
		my $tk = "$eleID\t$goid"; 
		defined $pair_ele_go{$tk} and next; 
		$pair_ele_go{$tk} = 1; 
		# Element ID to group ID; 
		push(@{$grpBg{'ele2grpID'}{$eleID}}, $goid); 
	}
	close($fh); 
	$grpBg{'eleNum'} = scalar(keys %{$grpBg{'eleID'}}); 
	for my $eleID (keys %{$grpBg{'ele2grpID'}}) {
		for my $goid (@{$grpBg{'ele2grpID'}{$eleID}}) {
			push(@{$grpBg{'grpID2ele'}{$goid}}, $eleID); 
		}
	}
	for my $goid (keys %{$grpBg{'grpID2ele'}}) {
		$grpBg{'grpSize_inGenome'}{$goid} = scalar( @{$grpBg{'grpID2ele'}{$goid}} ); 
	}
	return; 
}# load_goBg

sub load_subset {
	# Load data into %grpSub; Since the input file -inSubList has the same format, this function is general. 
	#   $grpSub{'ele2grpID'} = { eleID=>[grpID1, grpID2, ...], ... } 
	#   $grpSub{'grpID2ele'} = { grpID=>[eleID1, eleID2, ...], ... }
	#   $grpSub{'eleID'}     = { eleID=>1 }
	#   $grpSub{'eleNum'}    = scalar(keys %{$grpBg{'ele2grpID'}}) 
	#   $grpSub{'grpSize_inSubset'} = { grpID1=>NumberOfEleInGrpInSubset, grpID2=> , ... }
	#   
	my $fh = &openFH($gg{'inSubList'}, '<'); 
	while (<$fh>) {
		$. <= $gg{'skipHlineN'} and next; 
		chomp; 
		my @ta = &splitL("\t", $_); 
		( defined $ta[$gg{'colN'}] and $ta[$gg{'colN'}] !~ m!^\s*(NA|)\s*$!i ) or next; 
		my $eleID = $ta[$gg{'colN'}]; 
		defined $grpSub{'eleID'}{$eleID} and next; 
		$grpSub{'eleID'}{$eleID} = 1; 
		if (defined $grpBg{'ele2grpID'}{$eleID}) {
			for my $grpID (@{$grpBg{'ele2grpID'}{$eleID}}) {
				$grpSub{'grpSize_inSubset'}{$grpID} ++; 
				push(@{$grpSub{'ele2grpID'}{$eleID}}, $grpID); 
				push(@{$grpSub{'grpID2ele'}{$grpID}}, $eleID); 
			}
		} else {
			# This element has no group information. 
			# &stopErr("[Err] Failed to find element [$eleID] in background data.\n"); 
		}
	}
	close($fh); 
	$grpSub{'eleNum'} = scalar(keys %{$grpSub{'eleID'}}); 
	if (scalar(keys %{$grpSub{'grpSize_inSubset'}}) == 0) {
		&tsmsg("[Err] There is no $gg{'enrichType'} elements existing in subset!\n"); 
		exit(); 
	}
	for my $grpID (keys %{$grpSub{'grpID2ele'}}) {
		$grpSub{'grpSize_inSubset'}{$grpID} == scalar( @{$grpSub{'grpID2ele'}{$grpID}} ) or &stopErr("[Err] Different sizes for [$grpID].\n"); 
	}
	return; 
}# load_subset

sub run_fisherET {
	# Run fisher's exact test. 
	# Input should be %grpBg and %grpSub; 
	my $wrkDir = &fileSunhh::new_tmp_dir('create'=>1); 
	my $fh_oTab = &openFH("$wrkDir/input", '>'); 
	print {$fh_oTab} join("\t", qw/grpID n11 n12 n21 n22/)."\n"; 
	for my $grpID (keys %{$grpSub{'grpSize_inSubset'}}) {
		# In R - fisher.test() : 
		### total ele in genome : $grpBg{'eleNum'}                     $goBg{'eleNum'}
		### group ele in genome : $grpBg{'grpSize_inGenome'}{$grpID}   $goBg{'goid_cntInGenome'}{$goid}
		### total ele in subset : $grpSub{'eleNum'}                    $subInfo{'ele_cntInSubset'}
		### group ele in subset : $grpSub{'grpSize_inSubset'}{$grpID}  $subInfo{'goid_cntInSubset'}{$goid}
		my $n11 = $grpSub{'grpSize_inSubset'}{$grpID} ;  # in  subset & in  group(GO_ID); 
		my $n12 = $grpSub{'eleNum'} - $n11;              # in  subset & not group(GO_ID); 
		my $n21 = $grpBg{'grpSize_inGenome'}{$grpID} - $grpSub{'grpSize_inSubset'}{$grpID}; 
		                                                 # not subset & in  group(GO_ID); 
		my $n22 = $grpBg{'eleNum'} - $n11 - $n12 - $n21; # not subset & not group(GO_ID); 
		# $n11 + $n12 = ele_N_inSubset; 
		# $n21 + $n22 = ele_N_inRestset; 
		# $n11 + $n21 = ele_N_inGroup; 
		# $n12 + $n22 = ele_N_inRestgroup; 
		# $n11 + $n12 + $n21 + $n22 = ele_N_inGenomeBackground; 

		print {$fh_oTab} join("\t", "G$grpID", $n11, $n12, $n21, $n22)."\n"; 
	}
	close($fh_oTab); 
	my $fh_r = &openFH("$wrkDir/c.r", '>'); 
	print {$fh_r} <<"R1"; 
library(stats)
tbl     <- read.table( '$wrkDir/input', header=TRUE, sep="\\t", stringsAsFactors=FALSE, row.names=1 )
rawP    <- apply( tbl, MARGIN=1, FUN= function(x) { mat<-matrix( x[1:4], nrow=2, byrow=TRUE ); fisher.test(mat, alternative="greater")\$p.value } )
adjP    <- p.adjust( rawP, method= 'BH' ) # BH == FDR; 
ooo     <- data.frame( grpID=rownames(tbl), rawP=rawP, FDR=adjP )
ooo     <- ooo[ sort(ooo\$FDR, decreasing=T, index.return=TRUE)\$ix, ]
write.table( ooo, file='$wrkDir/output', append=FALSE, quote=FALSE, sep="\\t", row.names=FALSE, col.names=FALSE )

R1
	&exeCmd_1cmd("Rscript $wrkDir/c.r"); 

	my $fh_i = &openFH("$wrkDir/output", '<'); 
	my @out_lines; 
	if ( $gg{'enrichType'} eq 'go' ) {
		# For GO enrichment; 
		while (<$fh_i>) {
			chomp; 
			m!^\s*$! and next; 
			my @ta = &splitL("\t", $_); 
			my ($grpID, $rawP, $fdr) = @ta; 
			$grpID =~ s!^G!!; 
			$fdr < $gg{'padj'} or next; 
			push(@out_lines, 
			[
				$grpID, $fdr, $rawP, 
				$oboInfo{'goid2type'}{$grpID}, 
				$oboInfo{'goid2name'}{$grpID}, 
				$oboInfo{'goid2def'}{$grpID}, 
				join(";", @{$grpSub{'grpID2ele'}{$grpID}}), 
				join("_", 
					$grpSub{'grpSize_inSubset'}{$grpID}, 
					$grpSub{'eleNum'}, 
					$grpBg{'grpSize_inGenome'}{$grpID}, 
					$grpBg{'eleNum'}
				)
			]); 
		}
	} elsif ( $gg{'enrichType'} eq 'kegg' ) {
		# For KEGG enrichment; 
		while (<$fh_i>) {
			chomp; 
			m!^\s*$! and next; 
			my @ta = &splitL("\t", $_); 
			my ($grpID, $rawP, $fdr) = @ta; 
			$grpID =~ s!^G!!; 
			$fdr < $gg{'padj'} or next; 
			push(@out_lines, 
			[
				$grpID, $fdr, $rawP, 
				$grpBg{'grpID2def'}{$grpID}, 
				join(";", @{$grpSub{'grpID2ele'}{$grpID}}),
				join("_", 
					$grpSub{'grpSize_inSubset'}{$grpID}, 
					$grpSub{'eleNum'}, 
					$grpBg{'grpSize_inGenome'}{$grpID}, 
					$grpBg{'eleNum'}
				)
			]); 
		}
	}
	close($fh_i); 
	my $fh_o = &openFH("$gg{'opref'}.tab", '>'); 
	if ( $gg{'enrichType'} eq 'go' ) {
		# @out_lines > 0 and print {$fh_o} join("\t", qw/GO_ID FDR rawP GO_Type GO_name GO_def GO_eleID_sub sub_goN_totalN_bg_goN_totalN/)."\n"; 
		@out_lines > 0 and print {$fh_o} join("\t", qw/GO_ID FDR rawP GO_Type GO_name GO_def GO_eleID_sub subsetRatio backgroundRatio/)."\n"; 
		for (sort { $a->[3] cmp $b->[3] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @out_lines) {
			$_->[7] =~ m!^(\d+)_(\d+)_(\d+)_(\d+)$! or &stopErr("[Err] Bad format of sub_goN_totalN_bg_goN_totalN [$_->[7]]\n"); 
			my ($goN_sub, $totalN_sub, $goN_bg, $totalN_bg) = ($1,$2,$3,$4); 
			print {$fh_o} join("\t", @{$_}[0..6], "${goN_sub}_$totalN_sub", "${goN_bg}_$totalN_bg")."\n"; 
		}
	} elsif ( $gg{'enrichType'} eq 'kegg' ) {
		@out_lines > 0 and print {$fh_o} join("\t", qw/kegg_ID FDR rawP kegg_def eleID_sub subsetRatio backgroundRatio/)."\n"; 
		for (sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] || $a->[0] cmp $b->[0]} @out_lines) {
			$_->[5] =~ m!^(\d+)_(\d+)_(\d+)_(\d+)$! or &stopErr("[Err] Bad format of sub_goN_totalN_bg_goN_totalN [$_->[5]]\n"); 
			my ($goN_sub, $totalN_sub, $goN_bg, $totalN_bg) = ($1,$2,$3,$4); 
			print {$fh_o} join("\t", @{$_}[0..4], "${goN_sub}_$totalN_sub", "${goN_bg}_$totalN_bg")."\n"; 
		}
	}
	close($fh_o); 
	&fileSunhh::_rmtree($wrkDir); 
	return; 
}# run_fisherET() 

