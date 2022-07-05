#!/usr/bin/perl
use strict; 
use warnings; 
use File::Which; 
use LogInforSunhh; 
use fileSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"cpuN:i", 
	"path_fqC:s", 
	"path_montage:s", 
	"para_fqC:s", 
	"inFq:s@", 
	"inFqLis:s", 
	"outDir:s", 
	"noSummary!", 
	"onlySummary!", "pref_summ:s", 
); 
sub usage {
	print <<HH; 
################################################################################
# perl $0 -inFq in1.fq [-inFq in2.fq ...] -outDir dir_out
#
# -help 
# 
# -inFqLis   [] list of fq file names. 
# -path_fqC  [/path/to/fastqc]
# -para_fqC  [--nogroup --format fastq --kmers 5 --extract]
# 
# -cpuN      [10]
#
# -noSummary [] Generate summary if not given.
# -onlySummary   []
#
# -path_montage  [/path/to/montage] For summary usage. 
# -pref_summ     [pref]
#
################################################################################
HH
	exit 1; 
}

$opts{'help'} and &usage(); 
$opts{'path_fqC'} = $opts{'path_fqC'} // File::Which::which("fastqc"); 
-e $opts{'path_fqC'} or &stopErr("[Err] No fastqc file found at [$opts{'path_fqC'}].\n"); 
$opts{'path_montage'} = $opts{'path_montage'} // File::Which::which("montage"); 
( ! $opts{'noSummary'} ) and ( -e $opts{'path_montage'} or &stopErr("[Err] No montage exe file found at [$opts{'path_montage'}].\n") ); 
$opts{'pref_summ'} = $opts{'pref_summ'} // 'pref'; 
$opts{'para_fqC'} = $opts{'para_fqC'} // "--nogroup --format fastq --kmers 5 --extract"; 
$opts{'cpuN'} = $opts{'cpuN'} // 10; 
$opts{'outDir'} = $opts{'outDir'} // "dir_out"; 

my @inFqFiles; 
defined $opts{'inFq'} and push(@inFqFiles, @{$opts{'inFq'}}); 
if (defined $opts{'inFqLis'}) {
	my $fh = &openFH($opts{'inFqLis'}, '<'); 
	while (<$fh>) {
		m/^\s*#/ and next; 
		m/^\s*$/ and next; 
		s/\s+$//; 
		my @ta = split(/\s+/, $_); 
		$ta[0] =~ m/^\s*$/ and next; 
		push(@inFqFiles, $ta[0]); 
	}
	close($fh); 
}
scalar(@inFqFiles) > 0 or &stopErr("[Err] No input fastq files?\n"); 
my $inFq_string = join(" ", @inFqFiles); 

&tsmsg("[Rec] Total ", scalar(@inFqFiles), " input files to be processed.\n"); 

-d $opts{'outDir'} or mkdir($opts{'outDir'}, 0755); 
$opts{'onlySummary'} or &exeCmd("$opts{path_fqC} $opts{para_fqC} --threads $opts{'cpuN'} -o $opts{'outDir'} $inFq_string"); 

unless ( $opts{'noSummary'} ) {
	chdir($opts{'outDir'}); 
	my $sumDir = "00.Summary"; 
	-d $sumDir or mkdir($sumDir, 0755); 
	&exeCmd("$opts{path_montage} -label '\%d\%f' */*/per_base_quality.png             -tile 2x -frame 5 -shadow -geometry '+2+2>' $sumDir/$opts{pref_summ}_FqC_perBpQual.jpg"); 
	&exeCmd("$opts{path_montage} -label '\%d\%f' */*/sequence_length_distribution.png -tile 2x -frame 5 -shadow -geometry '+2+2>' $sumDir/$opts{pref_summ}_FqC_RdLenDistr.jpg"); 
	&exeCmd("$opts{path_montage} -label '\%d\%f' */*/per_sequence_gc_content.png      -tile 2x -frame 5 -shadow -geometry '+2+2>' $sumDir/$opts{pref_summ}_FqC_perSeqGC.jpg"); 
	&exeCmd("$opts{path_montage} -label '\%d\%f' */*/per_base_sequence_content.png    -tile 2x -frame 5 -shadow -geometry '+2+2>' $sumDir/$opts{pref_summ}_FqC_perBpSeqCont.jpg"); 
	&exeCmd("$opts{path_montage} -label '\%d\%f' */*/per_base_n_content.png           -tile 2x -frame 5 -shadow -geometry '+2+2>' $sumDir/$opts{pref_summ}_FqC_perBpNCont.jpg"); 
	&exeCmd("$opts{path_montage} -label '\%d\%f' */*/duplication_levels.png           -tile 2x -frame 5 -shadow -geometry '+2+2>' $sumDir/$opts{pref_summ}_FqC_DuplicLvl.jpg"); 
	&exeCmd("$opts{path_montage} -label '\%d\%f' */*/per_sequence_quality.png         -tile 2x -frame 5 -shadow -geometry '+2+2>' $sumDir/$opts{pref_summ}_FqC_perSeqQual.jpg"); 
}

&tsmsg("[Rec] All done.\n"); 


