#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use Parallel::ForkManager;
use File::Path qw(make_path remove_tree);
use File::Copy; 
use Cwd 'abs_path'; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"path_repClass:s", "para_RC:s", 
	"inFa:s", "tagUnknown:s", 
	"noChopID!", 
	"cpuN:i", 
	"out:s", 
); 

sub usage {
	print STDOUT <<HH; 
################################################################################
# perl $0 -inFa toBeClassified.fas -tagUnknown LTR
#
# -help
# -path_repClass [RepeatClassifier]
# -para_RC       [-engine ncbi]
#
# -out           [inFa.classified]
#
# -noChopID      [boolean] Do not chop ID if given. 
#
# -cpuN          [1] Parallelly run if larger than 1. 
# -pl_dealFa     [deal_fasta.pl] Must be given if cpuN > 1
#
################################################################################
HH
	exit 1; 
}
$opts{'path_repClass'} = $opts{'path_repClass'} // '/share/app/Annotation/repeatmodeler/RepeatModeler/RepeatClassifier'; 
$opts{'para_RC'} = $opts{'para_RC'} // '-engine ncbi'; 
$opts{'tagUnknown'} = $opts{'tagUnknown'} // 'LTR'; 
$opts{'cpuN'} = $opts{'cpuN'} // 1; 
$opts{'pl_dealFa'} = $opts{'pl_dealFa'} // `echo \$HOME/tools/github/NGS_data_processing/deal_fasta.pl`; 
chomp( $opts{'pl_dealFa'} ); 

$opts{'help'} and &usage(); 
defined $opts{'inFa'} and -e $opts{'inFa'} or &usage(); 
$opts{'out'} = $opts{'out'} // "$opts{'inFa'}.classified"; 

&tsmsg("[Rec] Begin $0\n"); 

my $tmp_file = 'tmp_toClass.fa'; 
unlink($tmp_file); 

unless ( $opts{'noChopID'} ) {
	open IN,'<',"$opts{'inFa'}" or &stopErr("[Err] $!\n"); 
	open OO,'>',"$tmp_file" or &stopErr("[Err] $1\n"); 
	while (<IN>) {
		chomp; 
		if (m/^>/) {
			m/^>(\S+)/ or &stopErr("[Err] Bad format of $_\n"); 
			$_ = ">$1"; 
			s!\#.*$!!; 
		}
		print OO "$_\n"; 
	}
	close OO; 
	close IN; 
} else {
	&exeCmd("cat $opts{'inFa'} > $tmp_file"); 
}

if ( $opts{'cpuN'} <= 1 ) {
	&exeCmd("$opts{'path_repClass'} -consensi $tmp_file $opts{'para_RC'}"); 
} else {
	-e "cut1_cutted" and File::Path::remove_tree( 'cut1_cutted', {'safe'=>1, 'keep_root'=>0} );
	&exeCmd("perl $opts{'pl_dealFa'} $tmp_file -cut 1 -cut_prefix cut1 1>cut1.std 2>cut1.err"); 
	-e "run_cut1" and File::Path::remove_tree( 'run_cut1', {'safe'=>1, 'keep_root'=>0} );
	mkdir('run_cut1'); 
	opendir DD,"cut1_cutted/" or &stopErr("[Err] $!\n");
	my @inFas; 
	for my $fn (readdir(DD)) {
		$fn =~ m/^\./ and next; 
		$fn =~ m/^cut1_\d+.fasta$/ or next; 
		mkdir("run_cut1/dir_${fn}/", 0755); 
		File::Copy::copy("cut1_cutted/$fn", "run_cut1/dir_${fn}/$fn") or &stopErr("[Err] Failed to copy file $fn\n"); 
		push(@inFas, $fn); 
	}
	closedir(DD); 
	@inFas = sort { my ($n1) = ($a =~ m/_(\d+)\.fasta$/); my ($n2) = ($b =~ m/_(\d+)\.fasta$/); $n1<=>$n2; } @inFas; 
	
	my $cwd = abs_path( File::Path::getcwd() ); 
	
	my $MAX_PROCESSES = $opts{'cpuN'} ; # Sometimes $parm{'cpuN'} - 1 may be better.
	my $pm = new Parallel::ForkManager($MAX_PROCESSES); 
	for my $fn (@inFas) {
		my $pid = $pm->start and next; 
		chdir("run_cut1/dir_${fn}/"); 
		&exeCmd("$opts{'path_repClass'} -consensi $fn $opts{'para_RC'}"); 
		chdir("../../"); 
		$pm->finish; 
	}
	$pm->wait_all_children; 
	
	chdir($cwd); 
	-e "${tmp_file}.classified" and unlink("${tmp_file}.classified"); 
	for my $fn (@inFas) {
		-e "run_cut1/dir_${fn}/${fn}.classified" or do { &tsmsg("[Err] File run_cut1/dir_${fn}/${fn}.classified not found.\n"); next; }; 
		&exeCmd("cat run_cut1/dir_${fn}/${fn}.classified >> ${tmp_file}.classified"); 
	}
}

open FF,'<',"${tmp_file}.classified" or &stopErr("[Err] $!\n"); 
open BB,'>',"$opts{'out'}" or &stopErr("[Err] $!\n"); 
while (<FF>) {
	chomp; 
	if (m/^>/) {
		s/^>(\S+)// or &stopErr("[Err] Bad format of $_\n"); 
		my $h1 = $1; 
		my $tag1 = ''; 
		if ( $h1 =~ s!\#(\S+)$!! ) {
			$tag1 = $1; 
			$tag1 eq 'Unknown' and $tag1 = $opts{'tagUnknown'}; 
		} else {
			&tsmsg("[Wrn] No tag information for ID [$h1]\n"); 
			$tag1 = $opts{'tagUnknown'}; 
		}
		$_ = ">${h1}#${tag1}$_"; 
		$_ =~ s/\s+$//; 
	}
	print BB "$_\n"; 
}
close BB; 
close FF; 

&tsmsg("[Rec] Finish $0\n"); 
