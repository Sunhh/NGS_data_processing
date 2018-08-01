#!/usr/bin/perl
# 2018-07-18 Use sortmerna to align reads to rRNA databases. 
# I used three LSU watermelon sequences from SILVA_132_LSURef_tax_silva_trunc to test the similarity between reference rDNA seqeuence and watermelon assemblies. 
# CMD : deal_fasta.pl /Data/Sunhh/database/db_fasta/rrna_silva/R132/SILVA_132_LSURef_tax_silva_trunc.Eukaryota.fa -res 'Citrullus lanatus' > t1.fa
# CMD : bn6 -query t1.fa -db db/watermelon_97103_pbChrV1.fa -evalue 1e-5 | awk ' $4 >= 300 ' | deal_table.pl -col_sort 2 > t1.fa.bn6.slct
# The hit regions in watermelon assembly were aligned to NCBI-Nt database to find best match sequence. 
# Then I think 85% is a good threshold for idenitty%, and I want to add 85% to coverage% requirement too. 
# I'd like to add E-value <= '1e-5' threshold too. 
# AGCB01009766.39202.42576        salGrp085       74.003  1981    367     110     121     2095    5693    3855    0.0     669     3375    7025    minus
#   Best match : Helianthus annuus uncharacterized LOC110926720 (LOC110926720), partial mRNA
# AGCB01009766.39202.42576        salGrp071       78.210  2625    470     76      753     3375    7       2531    0.0     1585    3375    10458   plus
#   Best match : Helianthus annuus fibroin heavy chain-like (LOC110905210), mRNA
# AGCB01005631.756.4253   WM97pbV1_Chr02  85.995  1221    118     24      1288    2487    36310379        36311567        0.0     1258    3498    37915939        plus
#   Best match : Cucumis sativus 18S ribosomal RNA gene, partial sequence; and internal transcribed spacer 1, 5.8S ribosomal RNA gene, internal transcribed spacer 2, and 26S ribosomal RNA gene, complete sequence
# AGCB01005631.756.4253           WM97pbV1_Chr09  86.392  316     42      1       1092    1407    16062986        16063300        2.34e-92        344     3498    37727573        plus
#   Best match : PREDICTED: Cucurbita maxima Eukaryotic 28S ribosomal RNA (LOC111475340), rRNA
# AGCB01009766.39202.42576        WM97pbV1_Chr09  86.392  316     42      1       689     1004    16062986        16063300        2.26e-92        344     3375    37727573   plus
#   Best match : PREDICTED: Cucurbita maxima Eukaryotic 28S ribosomal RNA (LOC111475340), rRNA
# AGCB01009767.1.2340     WM97pbV1_Chr02  87.437  1584    147     20      1       1582    36310588        36312121        0.0     1775    2340    37915939        plus
#   Best match : Cucumis sativus 18S ribosomal RNA gene, partial sequence; and internal transcribed spacer 1, 5.8S ribosomal RNA gene, internal transcribed spacer 2, and 26S ribosomal RNA gene, complete sequence
# AGCB01005631.756.4253   WM97pbV1_Chr07  87.850  428     34      8       2062    2472    13071911        13072337        3.67e-135       486     3498    31939013        plus
#   Best match : Symphoricarpos occidentalis voucher Steele 1358 26S ribosomal RNA gene, partial sequence
# AGCB01009767.1.2340     WM97pbV1_Chr08  90.326  827     76      3       51      876     7738399 7739222 0.0     1081    2340    28201227        plus
#   Best match : Cucumis sativus 18S ribosomal RNA gene, partial sequence; and internal transcribed spacer 1, 5.8S ribosomal RNA gene, internal transcribed spacer 2, and 26S ribosomal RNA gene, complete sequence
# AGCB01005631.756.4253   WM97pbV1_Chr06  90.390  385     19      6       2107    2473    19467231        19466847        2.84e-136       490     3498    29507460        minus
#   Best match : Symphoricarpos occidentalis voucher Steele 1358 26S ribosomal RNA gene, partial sequence
# AGCB01005631.756.4253   WM97pbV1_Chr01  90.416  1273    97      11      1223    2473    16010593        16011862        0.0     1652    3498    36935898        plus
#   Best match : Cucumis sativus 18S ribosomal RNA gene, partial sequence; and internal transcribed spacer 1, 5.8S ribosomal RNA gene, internal transcribed spacer 2, and 26S ribosomal RNA gene, complete sequence
#
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"inFq:s@", # Single file means RS, two files means R1/R2, others error; 
	"outPref:s@", # Output file prefix. 
	"db_sortmerna:s@", # ref,index_db:ref,index_db
	"db_sortmernaList:s", 
	"para_sortmerna:s",# -a 10 -m 8000 -e 1e-5 --otu_map --best 1 --id 0.85 --coverage 0.85 --log -v --blast '1 cigar qcov qstrand'
	"exe_sortmerna:s", # sortmerna 
	"help!", 
); 

my %gg; 
&setGlob(); 

for my $fq_aref (@{$gg{'inFq'}}) {
	my ($opref, $fq1, $fq2) = @$fq_aref; 
	my ($ofh1, $ofh2); 
	if (defined $fq2 and $fq2 ne '') {
		# paired-end read file
		my $cleanFq1 = &run_sortmerna( $fq1, "$gg{'tmpDir'}/r1", "${opref}_R1" ); 
		my $cleanFq2 = &run_sortmerna( $fq2, "$gg{'tmpDir'}/r2", "${opref}_R2" ); 
		$ofh1 = &openFH($cleanFq1, '>'); 
		$ofh2 = &openFH($cleanFq2, '>'); 
		&extract_fq_by_otu(["$gg{'tmpDir'}/r1.fq", "$gg{'tmpDir'}/r2.fq"], ["$gg{'tmpDir'}/r1_o_otus.txt", "$gg{'tmpDir'}/r2_o_otus.txt"], [$ofh1, $ofh2]); 
		close($ofh1); 
		close($ofh2); 
	} else {
		# Single-end read file 
		my $cleanFq1 = &run_sortmerna( $fq1, "$gg{'tmpDir'}/r1", "${opref}_RS" ); 
		$ofh1 = &openFH($cleanFq1, '>'); 
		&extract_fq_by_otu(["$gg{'tmpDir'}/r1.fq"], ["$gg{'tmpDir'}/r1_o_otus.txt"], [$ofh1]); 
		close($ofh1); 
	}
	&exeCmd_1cmd("rm -rf $gg{'tmpDir'}/*"); 
}
&fileSunhh::_rmtree($gg{'tmpDir'}); 

sub extract_fq_by_otu {
	my ($fqAR, $otuAR, $ofhAR) = @_; 
	my %id2RM; 
	for my $otuFile (@$otuAR) {
		open F1,'<',"$otuFile" or &stopErr("[Err] Failed to open file [$otuFile]\n"); 
		while (<F1>) {
			chomp; 
			my @ta = split(/\t/, $_); 
			shift(@ta); 
			for my $tb (@ta) {
				$id2RM{$tb} = 1; 
			}
		}
		close F1; 
	}

	if ( @$fqAR == 1 ) {
		open F2,'<',"$fqAR->[0]" or &stopErr("[Err] Failed to open fqAR[0] [$fqAR->[0]]\n"); 
		my ($l1, $l2, $l3, $l4); 
		while ($l1 = <F2>) {
			$l2 = <F2>; $l3 = <F2>; $l4 = <F2>; 
			chomp($l1); 
			$l1 =~ m!^\@(\S+)! or &stopErr("[Err] Bad read ID [$l1]\n"); 
			defined $id2RM{$1} and next; 
			print {$ofhAR->[0]} "$l1\n$l2$l3$l4"; 
		}
		close F2; 
	} elsif ( @$fqAR == 2 ) {
		open F2,'<',"$fqAR->[0]" or &stopErr("[Err] Failed to open fqAR[0] [$fqAR->[0]]\n"); 
		open F3,'<',"$fqAR->[1]" or &stopErr("[Err] Failed to open fqAR[1] [$fqAR->[1]]\n"); 
		my ($l1, $l2, $l3, $l4); 
		my ($r1, $r2, $r3, $r4); 
		while ($l1 = <F2>) {
			$l2 = <F2>; $l3 = <F2>; $l4 = <F2>; 
			$r1 = <F3>; $r2 = <F3>; $r3 = <F3>; $r4 = <F3>; 
			chomp($l1); 
			$l1 =~ m!^\@(\S+)! or &stopErr("[Err] Bad read ID [$l1]\n"); 
			defined $id2RM{$1} and next; 
			chomp($r1); 
			$r1 =~ m!^\@(\S+)! or &stopErr("[Err] Bad read ID [$r1]\n"); 
			defined $id2RM{$1} and next; 
			print {$ofhAR->[0]} "$l1\n$l2$l3$l4"; 
			print {$ofhAR->[1]} "$r1\n$r2$r3$r4"; 
		}
		close F2; 
		close F3; 
	} else {
		&stopErr("[Err] Bad input of fqAR [@$fqAR] in extract_fq_by_otu()\n"); 
	}
	
}# extract_fq_by_otu
sub run_sortmerna {
	my ($inFq, $resultPref, $cleanFqPref) = @_; 
	my $back_cleanFqName; 
	if      ( $inFq =~ m!\.gz$! ) {
		&exeCmd_1cmd("gzip -cd $inFq > ${resultPref}.fq") and &stopErr("[Err] Failed to gunzip $inFq\n"); 
		$back_cleanFqName = "${cleanFqPref}.fq.gz"; 
	} elsif ( $inFq =~ m!\.bz2$! ) {
		&exeCmd_1cmd("bzip2 -cd $inFq > ${resultPref}.fq") and &stopErr("[Err] Failed to bunzip2 $inFq\n"); 
		$back_cleanFqName = "${cleanFqPref}.fq.bz2"; 
	} else {
		&exeCmd_1cmd("ln -s $inFq ${resultPref}.fq") and &stopErr("[Err] Failed to link $inFq\n"); 
		$back_cleanFqName = "${cleanFqPref}.fq"; 
	}
	my $cmd = ""; 
	$cmd .= "$gg{'exe_sortmerna'} "; 
	$cmd .= " --ref $gg{'db_sortmerna'} "; 
	$cmd .= " $gg{'para_sortmerna'} "; 
	$cmd .= " --aligned ${resultPref}_o "; 
	$cmd .= " --reads ${resultPref}.fq "; 
	# I want to record stdout and stderr of sortmerna in screening;
	&exeCmd_1cmd($cmd) and &stopErr("[Err] Failed to run sortmerna: $cmd\n"); 
	return($back_cleanFqName); 
}# run_sortmerna() 

sub setGlob {
	$gg{'para_sortmerna'} = '-a 10 -m 8000 -e 1e-5 --otu_map --best 1 --id 0.85 --coverage 0.85 --log -v --blast \'1 cigar qcov qstrand\''; 
	$gg{'exe_sortmerna'}  = 'sortmerna'; 
	$gg{'oriDir'} = &fileSunhh::_abs_path_4link("./"); 
	$gg{'help_txt'} = <<HH; 
################################################################################
#  perl $0   -inFq in_R1.fq,in_R2.fq   -outPref outPrefix   -db_sortmerna db_fasta,db_indexdb
#
#  -help             [Bool]
#
#  -db_sortmernaList [filename] Used to shorten multiple -db_sortmerna .
#                      Format : db_fasta_1,db_indexdb_1
#                               db_fasta_2,db_indexdb_2
#                               ...
#
#  -para_sortmerna   ['$gg{'para_sortmerna'}']
#  -exe_sortmerna    ['$gg{'exe_sortmerna'}']
################################################################################
HH
	defined $opts{'help'} and &LogInforSunhh::usage($gg{'help_txt'}); 
	for my $t1 (qw/inFq outPref/) {
		defined $opts{$t1} or &stopErr("[Err] Parameter -$t1 is required.\n"); 
	}
	( defined $opts{'db_sortmerna'} or defined $opts{'db_sortmernaList'} ) or &stopErr("[Err] At least one of -db_sortmerna and -db_sortmernaList is required.\n"); 
	for (my $i=0; $i<@{$opts{'inFq'}}; $i++) {
		my $t1 = $opts{'inFq'}[$i]; 
		( defined $opts{'outPref'}[$i] and $opts{'outPref'}[$i] ne '' ) or &stopErr("[Err] -outPref should have the same number of -inFq\n"); 
		my @ta = map { s!^\s+|\s+$!!g; $_; } split(/,/, $t1); 
		@ta == 1 or @ta == 2 or &stopErr("[Err] Need one or two files in -inFq option. Now: |$t1|\n"); 
		push(@{$gg{'inFq'}}, [ &fileSunhh::_abs_path_4link($opts{'outPref'}[$i]) , map { &fileSunhh::_abs_path_4link($_) } @ta]); 
	}
	if (defined $opts{'db_sortmernaList'}) {
		my @tt1 = grep { $_ ne '' } map { $_->[0] } &fileSunhh::load_tabFile($opts{'db_sortmernaList'}); 
		push(@{$opts{'db_sortmerna'}}, @tt1); 
	}
	my %used_db; 
	for my $t1 (@{$opts{'db_sortmerna'}}) {
		$t1 =~ m!^\s*([^,\s]+),([^,\s]+)\s*$! or &stopErr("[Err] Bad input of -db_sortmerna [$t1]\n"); 
		my ($p1, $p2) = ($1, $2); 
		defined $used_db{"$p1\t$p2"} and next; 
		$used_db{"$p1\t$p2"} = 1; 
		$p1 = &fileSunhh::_abs_path_4link($p1); 
		$p2 = &fileSunhh::_abs_path_4link($p2); 
		push(@{$gg{'db_sortmerna_arr'}}, "$p1,$p2"); 
	}
	$gg{'db_sortmerna'} = join(':', @{$gg{'db_sortmerna_arr'}}); 
	defined $opts{'para_sortmerna'} and $gg{'para_sortmerna'} = $opts{'para_sortmerna'}; 
	defined $opts{'exe_sortmerna'} and $gg{'exe_sortmerna'} = $opts{'exe_sortmerna'}; 
	$gg{'tmpDir'} = &fileSunhh::_abs_path_4link(&fileSunhh::new_tmp_dir('create' => 1)); 
	return; 
}# setGlob() 



