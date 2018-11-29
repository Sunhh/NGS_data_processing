#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use mathSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"rdCntFn:s", 
	"tpmFn:s", 
	"rpkmFn:s", 
	"compareList:s", 
	"outFDRFn:s", 
	"skipRepeatSample!", 
	"exe_Rscript:s", 
	"addRPKM!", # Add mean value only. 
		"genLen:s", # Required if given -addRPKM
		"cutFDR:f", 
		"cutFC:f", 
		"minExp:f", 
	"pl_cnt2TPM:s", # /home/Sunhh/tools/github/NGS_data_processing/rnaseq_tools/cnvt_cnt_to_TPM.pl 
	"help!", 
); 

my %gg; 
&setup_glob(); 
&load_rdCnt(); 
&load_compareList(); 

for (my $i=0; $i<@{$gg{'compList'}}; $i++) {
	if      ( $gg{'compList'}[$i][0] eq 'exactTest' ) {
		&run_exactTest(@{$gg{'compList'}[$i]}); 
	} elsif ( $gg{'compList'}[$i][0] eq 'DESeq2' ) {
		&run_DESeq2(@{$gg{'compList'}[$i]}); 
	} else {
		&tsmsg("[Wrn] Skip unknown test Type [$gg{'compList'}[$i][0]]\n"); 
	}
}

&out_fdr(); 

sub out_fdr {
	my $fh1 = &openFH($gg{'outFDRFn'}, '>'); 
	for (@{$gg{'fdr_out'}}) {
		print {$fh1} join("\t", @$_)."\n"; 
	}
	close($fh1); 
}#out_fdr() 

sub _cnt2RPKM {
	my ($wrkDir, $g1_repN, $g2_repN) = @_; 

	my %expData; 
	# Output 4 columns (RPKMmean_1 RPKMmean_2 FoldChange FDR) instead of 1 column (FDR) for each comparison. 
	my $cmd_t1 = "perl $gg{'pl_cnt2TPM'} -rdCntFn $wrkDir/rdCnt -eleInfoFn $gg{'genLen'} -cntRPKM $wrkDir/rpkm > $wrkDir/tpm"; 
	&exeCmd_1cmd($cmd_t1) and &stopErr("[Err] Failed at CMD: $cmd_t1\n"); 
	my $fh2 = &openFH("$wrkDir/rpkm",'<'); 
	while (<$fh2>) {
		chomp; 
		my @ta=split(/\t/, $_); 
		if ($. == 1) {
			$ta[0] eq 'eleID' or &stopErr("[Err] Bad line 0 [$ta[0]] : $_\n"); 
			for (my $i=0; $i<$g1_repN; $i++) {
				$ta[$i+1] eq "G1_Rep$i" or &stopErr("[Err] Bad line $i : $_\n"); 
			}
			for (my $i=0; $i+$g1_repN+1<@ta; $i++) {
				$ta[$i+$g1_repN+1] eq "G2_Rep$i" or &stopErr("[Err] Bad line $i+$g1_repN : $_\n"); 
			}
			next; 
		}
		my $exp1 = &mathSunhh::_mean(@ta[ 1 .. $g1_repN ]); 
		my $exp2 = &mathSunhh::_mean(@ta[ ($g1_repN+1) .. $#ta ]); 
		my $fc; 
		if ( $gg{'minExp'} > 0 ) {
			if ( $exp1 >= $gg{'minExp'} and $exp2 >= $gg{'minExp'} ) {
				$fc = $exp2 / $exp1 ; 
			} elsif ( $exp1 >= $gg{'minExp'} ) {
				$fc = $gg{'minExp'} / $exp1 ; 
			} elsif ( $exp2 >= $gg{'minExp'} ) {
				$fc = $exp2 / $gg{'minExp'} ; 
			} else {
				$fc = 1; 
			}
		} else {
			$fc   = ( $exp1 > 0 ) ? ($exp2/$exp1) : ($exp2/0.01) ; 
		}
		$expData{$ta[0]} = [ $exp1, $exp2, $fc ]; 
	}
	close($fh2); 

	return(%expData); 
}# _cnt2RPKM() 

sub run_DESeq2 {
	my ($testType, $g1, $g2, $cmnDisp) = @_; 
	my $g1_ID = $g1->[0]; my @g1_idx = @{$g1->[1]}; 
	my $g2_ID = $g2->[0]; my @g2_idx = @{$g2->[1]}; 
	my $testID = 'ds.' . join("_VS_", $g1_ID, $g2_ID); 
	my $g1_repN = scalar(@g1_idx); 
	my $g2_repN = scalar(@g2_idx); 
	my $wrkDir = &fileSunhh::new_tmp_dir('create'=>1); 
	# Generate new read count file as input. 
	my $fh1 = &openFH("$wrkDir/rdCnt", '>'); 
	print {$fh1} join("\t", 'eleID', (map { "G1_Rep$_" } (0 .. $#g1_idx)), (map { "G2_Rep$_" } (0 .. $#g2_idx)))."\n"; 
	my $n = 0; 
	for my $line0 (@{$gg{'rdCnt'}}) {
		$n ++; 
		$n == 1 and next; 
		print {$fh1} join("\t", @{$line0}[0, @g1_idx, @g2_idx])."\n"; 
	}
	close($fh1); 
	my %expData; 
	$opts{'addRPKM'} and %expData = &_cnt2RPKM( $wrkDir, $g1_repN, $g2_repN ); 
	# Generate command R script. 
	my $fh3 = &openFH("$wrkDir/DESeq2.R", '>'); 
	print {$fh3} <<"DESeq2Test"; 
library(DESeq2)
library(dplyr)

### Load and filter read counts. 
rdCntFn   <- '$wrkDir/rdCnt'
rdCntData <- read.table( rdCntFn, header=T, sep=\"\\t\", row.names=1, stringsAsFactors=F )

rdCntData <- as.matrix(rdCntData)
# ( condition <- factor( gsub( "_Rep\\\\d", "", colnames(rdCntData)) ) )
( condition <- factor( c(rep(1, $g1_repN), rep(2, $g2_repN)) ) )
( coldata   <- data.frame(row.names=colnames(rdCntData), condition)  )

y1.dds <- DESeqDataSetFromMatrix(countData=rdCntData, colData=coldata, design=~condition)

# Filter by counts
y2.dds <- y1.dds[ rowSums(counts(y1.dds)) > 1, ]

y2.dds <- DESeq(y2.dds)
y2.res <- results(y2.dds)
y2.res.df <- data.frame( eleID=rownames(y2.res), as.data.frame(y2.res) )
# y2.res.df <- as.data.frame( y2.res )

### Output FDR. 
outFDR.ID    <- tibble::as_tibble( list(eleID=rownames(rdCntData)) )
toAddTbl     <- dplyr::as_tibble( y2.res.df ) %>% dplyr::select( "padj" ) %>% dplyr::mutate( eleID=rownames(y2.res.df) )
outTbl       <- dplyr::left_join( x=outFDR.ID, y=toAddTbl, by= "eleID" )
colnames(outTbl) <- c("eleID", "FDR")
write.table( outTbl, file="$wrkDir/fdr.tab", sep="\\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( y2.res.df, file="$wrkDir/DESeq2_res.tab", sep="\\t", col.names=TRUE, row.names=FALSE, quote=FALSE )

DESeq2Test

	# Run Rscript and get result. 
	&exeCmd_1cmd("$gg{'exe_Rscript'} $wrkDir/DESeq2.R"); 
	my @fdrTbl = &fileSunhh::load_tabFile( "$wrkDir/fdr.tab", 1 ); 
	@fdrTbl == @{$gg{'fdr_out'}} or &stopErr("[Err] fdrTbl num $#fdrTbl+1 != fdr_out num $#{$gg{'fdr_out'}}+1\n"); 
	if ( $opts{'addRPKM'} ) {
		push(@{$gg{'fdr_out'}[0]}, "exp.$g1_ID", "exp.$g2_ID", "FC.${g2_ID}.by.${g1_ID}", $testID, "DEG.$testID"); 
		for (my $i=1; $i<@{fdrTbl}; $i++) {
			$gg{'fdr_out'}[$i][0] eq $fdrTbl[$i][0] or &stopErr("[Err] Different ID $gg{'fdr_out'}[$i][0] / $fdrTbl[$i][0]\n"); 
			my $tag = &_degTag($gg{'cutFDR'}, $gg{'cutFC'}, $fdrTbl[$i][1], $expData{$fdrTbl[$i][0]}[2]); 
			push(@{$gg{'fdr_out'}[$i]}, @{$expData{$fdrTbl[$i][0]}}[0,1,2], $fdrTbl[$i][1], $tag); 
		}
	} else {
		push(@{$gg{'fdr_out'}[0]}, $testID); 
		for (my $i=1; $i<@{fdrTbl}; $i++) {
			$gg{'fdr_out'}[$i][0] eq $fdrTbl[$i][0] or &stopErr("[Err] Different ID $gg{'fdr_out'}[$i][0] / $fdrTbl[$i][0]\n"); 
			push(@{$gg{'fdr_out'}[$i]}, $fdrTbl[$i][1]); 
		}
	}
	&fileSunhh::_rmtree($wrkDir); 
	return; 
}# run_DESeq2() 

sub _degTag {
	my ($cutFDR, $cutFC, $fdr, $fc) = @_; 
	if      ( $fdr =~ m!^NA$!i ) {
		return('NA'); 
	} elsif ( $fdr >= $cutFDR ) {
		return('N'); 
	} elsif ( $cutFC <= 1 ) {
		if ( $fc > 1 ) {
			return('H'); 
		} elsif ( $fc < 1 ) {
			return('L'); 
		} else {
			return('N_lowExp'); 
		}
	} elsif ( $fc < 1/$cutFC ) {
		return('L'); 
	} elsif ( $fc > $cutFC  ) {
		return('H'); 
	} elsif ( $fc < 1 ) {
		return('L_lowFC'); 
	} elsif ( $fc > 1 ) {
		return('H_lowFC'); 
	} elsif ( $fc == 1 ) {
		return('N_lowExp'); 
	} else {
		&stopErr("[Err] Why here ? ($cutFDR, $cutFC, $fdr, $fc)\n"); 
	}
}# _degTag() 

sub run_exactTest {
	my ($testType, $g1, $g2, $cmnDisp) = @_; 
	my $g1_ID = $g1->[0]; my @g1_idx = @{$g1->[1]}; 
	my $g2_ID = $g2->[0]; my @g2_idx = @{$g2->[1]}; 
	my $testID = 'et.' . join("_VS_", $g1_ID, $g2_ID); 
	my $wrkDir = &fileSunhh::new_tmp_dir('create'=>1); 
	my $g1_repN = scalar(@g1_idx); 
	my $g2_repN = scalar(@g2_idx); 
	# Generate new read count file as input. 
	my $fh1 = &openFH("$wrkDir/rdCnt", '>'); 
	print {$fh1} join("\t", 'eleID', (map { "G1_Rep$_" } (0 .. $#g1_idx)), (map { "G2_Rep$_" } (0 .. $#g2_idx)))."\n"; 
	my $n = 0; 
	for my $line0 (@{$gg{'rdCnt'}}) {
		$n ++; 
		$n == 1 and next; 
		print {$fh1} join("\t", @{$line0}[0, @g1_idx, @g2_idx])."\n"; 
	}
	close($fh1); 
	my %expData; 
	$opts{'addRPKM'} and %expData = &_cnt2RPKM( $wrkDir, $g1_repN, $g2_repN ); 

	# Generate command R script. 
	my $code_test = <<"ET_TEST"; 
y2.et.estD   <- estimateDisp( y2, design=y2.et.design )
y2.et.et     <- exactTest( y2.et.estD )
ET_TEST
	if (@g1_idx == 1 and @g2_idx == 1) {
		$code_test = <<"ET_TEST_NOREP"; 
y2.et.et     <- exactTest( y2, dispersion=$cmnDisp )
ET_TEST_NOREP
	}
	my $fh2 = &openFH("$wrkDir/edgeR_exactTest.R", '>'); 
	print {$fh2} <<"RexactTest"; 
library(dplyr)
library(edgeR)

### Load and filter read counts. 
rdCntFn   <- '$wrkDir/rdCnt'
rdCntData <- read.table( rdCntFn, header=T, sep=\"\\t\", row.names=1, stringsAsFactors=F )
y1.grp    <- factor( c(rep(1, $g1_repN), rep(2, $g2_repN)) )
y1        <- DGEList( rdCntData, group=y1.grp )
y1        <- calcNormFactors(y1)
y2.keep   <- rowSums(cpm(y1)>1) >= 2
y2        <- y1[ y2.keep, , keep.lib.sizes=FALSE ]
y2        <- calcNormFactors(y2)
y2.grp    <- y1.grp

### Compute FDR. 
y2.et.design <- model.matrix( ~y2.grp )
$code_test
y2.et.topTag <- topTags( y2.et.et, n=Inf )
# table( y2.et.topTag\$table\$FDR < 0.01 ) 

### Output FDR. 
outFDR.ID    <- tibble::as_tibble( list(eleID=rownames(rdCntData)) )
toAddTbl     <- tibble::as_tibble( y2.et.topTag\$table ) %>% dplyr::select( "FDR" ) %>% dplyr::mutate( eleID=rownames(y2.et.topTag\$table) )
outTbl       <- dplyr::left_join( x=outFDR.ID, y=toAddTbl, by= "eleID" )
write.table( outTbl, file="$wrkDir/fdr.tab", sep="\\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( y2.et.topTag\$table, file="$wrkDir/topTags.tab", sep="\\t", col.names=TRUE, row.names=TRUE, quote=FALSE )

RexactTest
	close($fh2); 
	# Run Rscript and get result. 
	&exeCmd_1cmd("$gg{'exe_Rscript'} $wrkDir/edgeR_exactTest.R"); 
	my @fdrTbl = &fileSunhh::load_tabFile( "$wrkDir/fdr.tab", 1 ); 
	@fdrTbl == @{$gg{'fdr_out'}} or &stopErr("[Err] fdrTbl num $#fdrTbl+1 != fdr_out num $#{$gg{'fdr_out'}}+1\n"); 

	if ( $opts{'addRPKM'} ) {
		push(@{$gg{'fdr_out'}[0]}, "exp.$g1_ID", "exp.$g2_ID", "FC.${g2_ID}.by.${g1_ID}", $testID, "DEG.$testID"); 
		for (my $i=1; $i<@{fdrTbl}; $i++) {
			$gg{'fdr_out'}[$i][0] eq $fdrTbl[$i][0] or &stopErr("[Err] Different ID $gg{'fdr_out'}[$i][0] / $fdrTbl[$i][0]\n"); 
			my $tag = &_degTag($gg{'cutFDR'}, $gg{'cutFC'}, $fdrTbl[$i][1], $expData{$fdrTbl[$i][0]}[2]); 
			push(@{$gg{'fdr_out'}[$i]}, @{$expData{$fdrTbl[$i][0]}}[0,1,2], $fdrTbl[$i][1], $tag); 
		}
	} else {
		push(@{$gg{'fdr_out'}[0]}, $testID); 
		for (my $i=1; $i<@{fdrTbl}; $i++) {
			$gg{'fdr_out'}[$i][0] eq $fdrTbl[$i][0] or &stopErr("[Err] Different ID $gg{'fdr_out'}[$i][0] / $fdrTbl[$i][0]\n"); 
			push(@{$gg{'fdr_out'}[$i]}, $fdrTbl[$i][1]); 
		}
	}
	&fileSunhh::_rmtree($wrkDir); 
	return; 
}# run_exactTest() 


sub load_compareList {
	# Load compare list to $gg{'compList'}
	my $fh = &openFH($gg{'compareList'}, '<'); 
	while (<$fh>) {
		chomp; 
		my @ta = &splitL("\t", $_); 
		if ($ta[0] =~ m!^\s*(exactTest|DESeq2)\s*$!i) {
			# Only two groups; 
			my (%grp1, %grp2); 
			for my $a1 (split(/;/, $ta[1])) {
				$a1 =~ s!^\s+|\s+$!!g; 
				$grp1{$a1} = 1; 
			}
			for my $a2 (split(/;/, $ta[2])) {
				$a2 =~ s!^\s+|\s+$!!g; 
				$grp2{$a2} = 1; 
			}
			# $ta[1] eq $ta[2] and &stopErr("[Err] Same groups $ta[1]\n"); 
			my (@i1, @i2); 
			for (my $i=0; $i<@{$gg{'grpInfo'}}; $i++) {
				if      ( defined $grp1{ $gg{'grpInfo'}[$i][0] } ) {
					push(@i1, $gg{'grpInfo'}[$i][2]); 
				} elsif ( defined $grp2{ $gg{'grpInfo'}[$i][0] } ) {
					push(@i2, $gg{'grpInfo'}[$i][2]); 
				}
			}
			my $repN_1 = scalar(@i1); 
			my $repN_2 = scalar(@i2); 
			$repN_1 == 0 and do { &tsmsg("[Wrn] No sample found for group [$ta[1]] in line : $_\n"); next; }; 
			$repN_2 == 0 and do { &tsmsg("[Wrn] No sample found for group [$ta[2]] in line : $_\n"); next; }; 
			my $cmnDisp = ''; 
			defined $ta[3] and $ta[3] !~ m!^(NA|\s*)$!i and $cmnDisp = $ta[3]; 
			if ( $repN_1 == 1 and $repN_2 == 1 and $cmnDisp eq '') {
				&tsmsg("[Wrn] I need a common dispersion if there are no replicates in any group. Skip line: $_\n"); 
				next; 
			}
			# ( @i1 >= 2 and @i2 >= 2 ) or &stopErr("[Err] Not enough replicates for $_.\n"); 
			&tsmsg("[Msg] Loading $ta[0] between group 1 [$ta[1]] x$repN_1 and group 2 [$ta[2]] x$repN_2\n"); 
			if ( $ta[0] =~ m!^\s*exactTest\s*$!i ) {
				push(@{$gg{'compList'}}, ['exactTest', [$ta[1], [@i1]], [$ta[2], [@i2]], $cmnDisp]); # ([testType, [grp1_ID, [grp1_idx]], [grp2_ID, [grp2_idx]], commnDisper], [testType, [], [], cD, cD])
			} elsif ( $ta[0] =~ m!^\s*DESeq2\s*$! ) {
				push(@{$gg{'compList'}}, ['DESeq2', [$ta[1], [@i1]], [$ta[2], [@i2]], $cmnDisp]); # ([testType, [grp1_ID, [grp1_idx]], [grp2_ID, [grp2_idx]], commnDisper], [testType, [], [], cD, cD])
			} else {
				&stopErr("[Err] Unknown [$ta[0]]\n"); 
			}
		} else {
			&tsmsg("[Wrn] Skip unparsable line : $_\n"); 
		}
	}
	close ($fh); 
	return; 
}# load_compareList() 

sub load_rdCnt {
	# Load rdCntFn to $gg{'rdCnt'}, $gg{'grpInfo'}
	$gg{'rdCnt'} = [ &fileSunhh::load_tabFile( $gg{'rdCntFn'}, 1 ) ]; 
	my %usedID; 
	for (my $i=1; $i<@{$gg{'rdCnt'}[0]}; $i++) {
		$gg{'rdCnt'}[0][$i] =~ m!^(\S+)_Rep(\d+)!i or &stopErr("[Err] Bad input sample_ID [$gg{'rdCnt'}[0][$i]] which should be [\\S+_Rep\\d+]\n"); 
		if ( defined $usedID{$gg{'rdCnt'}[0][$i]} ) {
			if ( $opts{'skipRepeatSample'} ) {
				next; 
			} else {
				&stopErr("[Err] There are same sample_IDs in rdCnt table [$gg{'rdCnt'}[0][$i]]\n"); 
			}
		}
		push(@{$gg{'grpInfo'}}, [$1, $2, $i, $gg{'rdCnt'}[0][$i]]); # ([grpID, repN, colN, sampleID], [], ...)
		$usedID{$gg{'rdCnt'}[0][$i]} = 1; 
	}
	@{$gg{'grpInfo'}} = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } @{$gg{'grpInfo'}}; 
	for (@{$gg{'rdCnt'}}) {
		push(@{$gg{'fdr_out'}}, [$_->[0]]); 
	}
	return; 
}# load_rdCnt 

sub setup_glob {
	$gg{'exe_Rscript'} //= 'Rscript'; 
	$gg{'pl_cnt2TPM'}  //= '/home/Sunhh/tools/github/NGS_data_processing/rnaseq_tools/cnvt_cnt_to_TPM.pl'; 
	$gg{'cutFDR'}      //= 0.01; 
	$gg{'cutFC'}       //= 2; 
	$gg{'minExp'}      //= 0.01; 
	$gg{'htxt'} = <<"H1"; 
################################################################################
# perl $0 -rdCntFn rdCnt.Both.use -compareList comp_grpPair_list -outFDRFn outFDR.list
#
# -help 
#
# -rdCntFn        [filename] required. Example : 
#                   EleID                   G1_Rep1 G1_Rep2 G1_Rep3 G2_Rep1 G2_Rep2
#                   Cla97C00G000010.1       497     819     550     842     1104
#                   Cla97C00G000020.1       510     684     581     463     742
#                   Cla97C00G000030.1       44      83      50      70      56
#                   Cla97C00G000040.1       433     645     408     487     469
#                   Cla97C00G000050.1       189     354     291     195     261
#
# -compareList     [filename] required. Example : 
#                   exactTest \\t G1 \\t G2 \\n
#                   exactTest \\t G1 \\t G3 \\n
#                   DESeq2    \\t G1 \\t G2 \\n
#
# -exe_Rscript     ['$gg{'exe_Rscript'}']; 
# -pl_cnt2TPM      [$gg{'pl_cnt2TPM'}]
# 
# -addRPKM         [Boolean]
#   -genLen        [filename] for -eleInfoFn in -pl_cnt2TPM ; 
#   -cutFDR
#   -cutFC
#
################################################################################
H1
	for my $k1 (qw/rdCntFn compareList outFDRFn/) {
		defined $opts{$k1} or &LogInforSunhh::usage($gg{'htxt'}); 
		$gg{$k1} = $opts{$k1}; 
	}
	for my $k2 (qw/pl_cnt2TPM exe_Rscript genLen cutFDR cutFC minExp/) {
		defined $opts{$k2} and $gg{$k2} = $opts{$k2}; 
	}
	for my $k3 (qw/exe_Rscript/) {
		$gg{$k3} = &fileSunhh::_abs_path_4link(&fileSunhh::_which($gg{$k3})); 
	}
}# setup_glob() 


