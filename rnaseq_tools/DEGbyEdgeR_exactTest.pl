#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
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
	"help!", 
); 

my %gg; 
&setup_glob(); 
&load_rdCnt(); 
&load_compareList(); 

for (my $i=0; $i<@{$gg{'compList'}}; $i++) {
	if ($gg{'compList'}[$i][0] eq 'exactTest') {
		&run_exactTest(@{$gg{'compList'}[$i]}); 
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

sub run_exactTest {
	my ($testType, $g1, $g2) = @_; 
	my $g1_ID = $g1->[0]; my @g1_idx = @{$g1->[1]}; 
	my $g2_ID = $g2->[0]; my @g2_idx = @{$g2->[1]}; 
	my $testID = 'et.' . join("_VS_", $g1_ID, $g2_ID); 
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
	# Generate command R script. 
	my $fh2 = &openFH("$wrkDir/edgeR_exactTest.R", '>'); 
	print {$fh2} <<"RexactTest"; 
library(dplyr)
library(edgeR)

### Load and filter read counts. 
rdCntFn   <- '$wrkDir/rdCnt'
rdCntData <- read.table( rdCntFn, header=T, sep=\"\\t\", row.names=1, stringsAsFactors=F )
y1.grp    <- factor( gsub( "_Rep\\\\d", "", colnames(rdCntData)) )
y1        <- DGEList( rdCntData, group=y1.grp )
y1        <- calcNormFactors(y1)
y2.keep   <- rowSums(cpm(y1)>1) >= 2
y2        <- y1[ y2.keep, , keep.lib.sizes=FALSE ]
y2        <- calcNormFactors(y2)
y2.grp    <- y1.grp

### Compute FDR. 
y2.et.design <- model.matrix( ~y2.grp )
y2.et.estD   <- estimateDisp( y2, design=y2.et.design )
y2.et.et     <- exactTest( y2.et.estD )
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
	push(@{$gg{'fdr_out'}[0]}, $testID); 
	for (my $i=1; $i<@{fdrTbl}; $i++) {
		$gg{'fdr_out'}[$i][0] eq $fdrTbl[$i][0] or &stopErr("[Err] Different ID $gg{'fdr_out'}[$i][0] / $fdrTbl[$i][0]\n"); 
		push(@{$gg{'fdr_out'}[$i]}, $fdrTbl[$i][1]); 
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
		if ($ta[0] =~ m!^\s*exactTest\s*$!i) {
			# Only two groups; 
			$ta[1] =~ s!^\s+|\s+$!!g; 
			$ta[2] =~ s!^\s+|\s+$!!g; 
			$ta[1] eq $ta[2] and &stopErr("[Err] Same groups $ta[1]\n"); 
			my (@i1, @i2); 
			for (my $i=0; $i<@{$gg{'grpInfo'}}; $i++) {
				if (     $gg{'grpInfo'}[$i][0] eq $ta[1]) {
					push(@i1, $gg{'grpInfo'}[$i][2]); 
				} elsif ($gg{'grpInfo'}[$i][0] eq $ta[2]) {
					push(@i2, $gg{'grpInfo'}[$i][2]); 
				}
			}
			( @i1 >= 2 and @i2 >= 2 ) or &stopErr("[Err] Not enough replicates for $_.\n"); 
			&tsmsg("[Msg] Loading exactTest between group 1 [$ta[1]] and group 2 [$ta[2]]\n"); 
			push(@{$gg{'compList'}}, ['exactTest', [$ta[1], [@i1]], [$ta[2], [@i2]]]); # ([testType, [grp1_ID, [grp1_idx]], [grp2_ID, [grp2_idx]] ], [testType, [], []])
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
#
# -exe_Rscript     ['$gg{'exe_Rscript'}']; 
#
################################################################################
H1
	for my $k1 (qw/rdCntFn compareList outFDRFn/) {
		defined $opts{$k1} or &LogInforSunhh::usage($gg{'htxt'}); 
		$gg{$k1} = $opts{$k1}; 
	}
	defined $opts{'exe_Rscript'} and $gg{'exe_Rscript'} = $opts{'exe_Rscript'}; 
	$gg{'exe_Rscript'} = &fileSunhh::_which($gg{'exe_Rscript'}); 
}# setup_glob() 


