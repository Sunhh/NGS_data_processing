#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use SeqAlnSunhh; 
my $sas = SeqAlnSunhh->new(); 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"db:s", # /data/Sunhh/database/db_bwa/Watermelon/WM97_v6.chr.fa
	"exe_tophat:s", # tophat
	"outDir:s", # ./tophat_out
	"para_tophat:s", # -p $cpuN --library-type=fr-firststrand --read-mismatches 1 --splice-mismatches 0 --min-intron-length 30
	"printCmd!", 

	"inFq1:s", "inFq2:s", "alnTogether!", 
	"para_list:s", 
	"help!", 
); 

sub usage {
	print <<HH;
################################################################################
#  perl $0 
#
# -para_list  [''] A file recording all information needed in bwa alignment. 
#                  This will overwrite the general settings. 
#
#  -inFq1     input_R1 fastq
#  -inFq2     input_R2 fastq
#  -db        bwa database indexed
#  -outDir    [./tophat_out]
#
#  -exe_tophat         ['tophat']
#
#  -para_tophat        [''] used in 'tophat'. Could be '-p \$cpuN --library-type=fr-firststrand --read-mismatches 1 --splice-mismatches 0 --min-intron-length 30'
#
#  -alnTogether        [FALSE] Run tophat only one time with all input read files together. 
#  
#  -printCmd        [FALSE] Only output commands if given. 
################################################################################
HH
	exit; 
}

$opts{'help'} and &usage(); 

for (qw/db/) {
	$opts{$_} //= 'NA'; 
}
$opts{'index'} = $opts{'db'}; 

$opts{'exe_tophat'}  //= 'tophat'; 
$opts{'printCmd'}    //= 0; 
$opts{'outDir'}      //= './tophat_out'; 
$opts{'para_tophat'} //= ''; 

my @batch; 
if (defined $opts{'inFq1'}) {
	push(@batch, { %opts }); 
}
if ( defined $opts{'para_list'} ) {
	my @header; 
	open F,'<',"$opts{'para_list'}" or die; 
	my $hh = <F>; 
	chomp($hh); 
	@header = split(/\t/, $hh); 
	while (<F>) {
		chomp; m/^(\s*$|\s*#)/ and next; 
		my @ta = split(/\t/, $_); 
		my %cur_opts; 
		for ( my $i=0; $i<@header; $i++ ) {
			$cur_opts{$header[$i]} = $ta[$i]; 
		}
		push(@batch, { %cur_opts }); 
		for my $tk ( keys %opts ) {
			defined $batch[-1]{$tk} or $batch[-1]{$tk} = $opts{$tk}; 
		}
	}
	close F; 
}

$#batch == -1 and &usage(); 
if ($opts{'alnTogether'}) {
	my (@inFqPE1, @inFqPE2); 
	for ( my $i=0; $i<@batch; $i++ ) {
		$batch[$i]{'inFq1'} //= ''; 
		$batch[$i]{'inFq2'} //= ''; 
		push(@inFqPE1, $batch[$i]{'inFq1'}); 
		push(@inFqPE2, $batch[$i]{'inFq2'}); 
	}
	my %parm; 
	for ( keys %{$batch[0]} ) {
		$_ =~ m/^inFq[12]$/ and next; 
		$parm{$_} = $batch[0]{$_}; 
	}
	$sas->aln_tophat2( 'inFq1'=>\@inFqPE1, 'inFq2'=>\@inFqPE2, %parm ); 
} else {
	for (my $i=0; $i<@batch; $i++) {
		$batch[$i]{'inFq1'} //= ''; 
		$batch[$i]{'inFq2'} //= ''; 
		$batch[$i]{'inFq1'} = [ map { s!^\s+|\s+$!!g; $_ } split(/,/, $batch[$i]{'inFq1'}) ]; 
		$batch[$i]{'inFq2'} = [ map { s!^\s+|\s+$!!g; $_ } split(/,/, $batch[$i]{'inFq2'}) ]; 
		$sas->aln_tophat2( %{$batch[$i]} ); 
		# $sas->aln_tophat2( 'inFq1'=>[ $batch[$i]{'inFq1'} ], 'inFq2'=>[ $batch[$i]{'inFq2'} ], 'index'=>$batch[$i]{'db'}, 'outDir'=>$batch[$i]{'outDir'}, 'para_tophat'=>$batch[$i]{'para_tophat'}, 'printCmd'=>$batch[$i]{'printCmd'}); 
	}
}


&tsmsg("[Rec]All done.\n"); 
