#!/usr/bin/perl
use strict; 
use warnings; 
use Cwd 'abs_path'; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"aln_type:s", # prot2genome
	"para_spaln:s", "cpuN:i", 
	"inFa:s@", "inFaLis:s", 
	"db:s", "needIndex!", 
	"cnvt2maker!", "pl_cnvt2maker:s", 
	"outMerge:s", 
	"printCmd!", 
	"help!", 
); 

sub usage {
	print <<HH;
################################################################################
# perl $0 -inFa protein.fa [-inFa in_prot_2.fa] -db db_genome
# 
# -out             [inFa.spaln.gff3]
# -outMerge        [''] If given, all gff3 outputs will be merged and treated with 'gff3_merge' script. 
#                       This will force trigger -cnvt2maker
# -inFaLis         ['']
# -aln_type        ['prot2genome'] Could be 'est2genome'
# -para_spaln      [Different according to aln_type]
# -cpuN            [1] 
# 
# -needIndex       [FALSE] The tmp-index files will be removed after used. 
# -printCmd        [FALSE]
# 
# -cnvt2maker      [FALSE] Convert output file to maker format if given. 
# -pl_cnvt2maker   [\$HOME/tools/github/NGS_data_processing/annot_tools/cnvt_spaln2makerAln_prot_gff3.pl]
################################################################################
HH
	exit 1; 
}

$opts{'help'} and &usage(); 
$opts{'printCmd'} //= 0; 
my $pl_cnvt2maker = $opts{'pl_cnvt2maker'} // `echo \$HOME/tools/github/NGS_data_processing/annot_tools/cnvt_spaln2makerAln_prot_gff3.pl`; 
chomp($pl_cnvt2maker); 
defined $opts{'outMerge'} and $opts{'cnvt2maker'} = 1; 

# Record input files
if (defined $opts{'inFaLis'}) {
	open F,'<',"$opts{'inFaLis'}" or die; 
	while (<F>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		push(@{$opts{'inFa'}}, $ta[0]); 
	}
	close F; 
}
$opts{'inFa'} // &usage(); 
$opts{'db'}   // &usage(); 
my $dbID = $opts{'db'}; 

# If -needIndex 
my $ori_dir = abs_path('.'); 
my $aln_dbs = `echo \$ALN_DBS`; 
chomp($aln_dbs); 
if ( $opts{'needIndex'} ) {
	my $db_path = abs_path($opts{'db'}); 
	$aln_dbs eq '' and &stopErr("[Err] No \$ALN_DBS assigned in environment.\n"); 
	my $tmp_dbID = 'tmpDB0001'; 
	while ( -e "$tmp_dbID.mfa" ) {
		$tmp_dbID ++; 
	}
	&exeCmd_1cmd("cd $aln_dbs; ln -s $db_path ./$tmp_dbID.mfa; ./makeidx.pl -inp $tmp_dbID.mfa; cd -;", $opts{'printCmd'}); 
	&tsmsg("[Rec]Change database ID from [$dbID] to [$tmp_dbID]\n"); 
	$dbID = $tmp_dbID; 
}

# Set default align parameters 
$opts{'aln_type'} //= 'prot2genome'; 
$opts{'aln_type'} = lc( $opts{'aln_type'} ); 
if ( $opts{'aln_type'} eq 'prot2genome' ) {
	$opts{'para_spaln'} //= ' -t1 -M4 -Q7 -O0 -LS -ya2'; 
} elsif ( $opts{'aln_type'} eq 'est2genome' ) {
	$opts{'para_spaln'} //= ' -t1 -M4 -Q7 -O0 -LS -ya2'; 
} else {
	&stopErr("[Err] Unknown -aln_type [$opts{'aln_type'}]\n"); 
}
if ( defined $opts{'cpuN'}  ) {
	$opts{'para_spaln'} =~ s!(^|\s)\-t\d+(\s|$)!$1$2!; 
	$opts{'para_spaln'} = "-t$opts{'cpuN'} $opts{'para_spaln'}"; 
}

# Run spaln
my @oFiles; 
for my $inFa ( @{$opts{'inFa'}} ) {
	$inFa = abs_path($inFa); 
	my $outFile = $opts{'out'} // "$inFa.spaln.gff3";
	&exeCmd_1cmd("spaln $opts{'para_spaln'} -o $outFile -d$dbID $inFa", $opts{'printCmd'}); 
	if ($opts{'cnvt2maker'}) {
		if ( $opts{'aln_type'} eq 'prot2genome') {
			&exeCmd_1cmd("perl $pl_cnvt2maker $outFile > $outFile.maker", $opts{'printCmd'}); 
		} elsif ( $opts{'aln_type'} eq 'est2genome' ) {
			&exeCmd_1cmd("perl $pl_cnvt2maker $outFile | perl -pe ' s!\\tprotein_match\\t!\\tmatch\\t!; ' > $outFile.maker", $opts{'printCmd'}); 
		} else {
			&exeCmd_1cmd("perl $pl_cnvt2maker $outFile > $outFile.maker", $opts{'printCmd'}); 
		}
		&exeCmd_1cmd("mv $outFile.maker $outFile", $opts{'printCmd'}); 
	}
	push(@oFiles, $outFile); 
}
if ( defined $opts{'outMerge'} ) {
	my $in_str = join(' ', @oFiles); 
	&exeCmd_1cmd("gff3_merge -l -o $opts{'outMerge'} $in_str", $opts{'printCmd'}); 
}
if ( $opts{'needIndex'} ) {
	my @idx_files; 
	for my $suff (qw/.mfa .seq .idx .grp .ent .bkn .bkp/) {
		push(@idx_files, "${dbID}$suff"); 
	}
	my $idx_str = join(' ', @idx_files); 
	&exeCmd_1cmd("cd $aln_dbs ; rm $idx_str ; cd -", $opts{'printCmd'}); 
}
