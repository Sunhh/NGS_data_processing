#!/usr/bin/perl
# 2018-10-26 : step X2 . This is used to remove duplicated reads in pref_comb.nonTrans.bam, which is output of stepX2; 
#   Input  : -pref , -wrk_dir , -inBam , -exe_samtools
#   Output : pref_comb.nonTransRmdup.bam
use strict; 
use warnings; 
use LogInforSunhh;
use SeqAlnSunhh;
use fileSunhh;
use Getopt::Long;
my %opts;
GetOptions(\%opts, 
	"help!", 
	"pref:s",    # Project prefix; 
	"wrk_dir:s", # './', Assign a directory to store all resulting bam files, instead of in current folder.
	"inBam:s",   # Input bam file, which should be $gg{'wrk_dir'}/$gg{'pref'}_comb.nonTrans.bam
	"outBam:s",  # The combined bam file output. This will replace the output file "$gg{'wrk_dir'}/$gg{'pref'}_comb.nonTransRmdup.bam" 

	"exe_samtools:s", 
); 

my %gg; 
&setGlob(); 
&applyOpt(); 
&stepX2_dedup(); 

sub stepX2_dedup {
	# Step X2 : Remove duplicated reads in pref_comb.nonTrans.bam; 
	# Produce : $gg{'wrk_dir'}/$gg{'pref'}_comb.nonTransRmdup.bam
	
	-e $gg{'inBam'} or &stopErr("[Err] -inBam $gg{'inBam'} not found. \n"); 
	open F,'-|', "$gg{'exe_samtools'} view -h $gg{'inBam'} " or &stopErr("[Err] Failed to read bam file [$gg{'inBam'}]\n"); 
	open O,'|-', "$gg{'exe_samtools'} view -o $gg{'outBam'} - " or &stopErr("[Err] Failed to write bam file [$gg{'outBam'}]\n"); 
	my %h; 
	while (<F>) {
		m!^\@! and do { print O $_; next; }; 
		chomp; 
		my @ta = split(/\t/, $_); 
		my $tk = join("\t", @ta[2,3,5,9]); 
		defined $h{$tk} and next; 
		$h{$tk} = 1; 
		print O "$_\n"; 
	}
	close O; 
	close F; 

	return; 
}# stepX2_dedup() 

sub applyOpt {
	$opts{'help'} and &LogInforSunhh::usage($gg{'help_txt'}); 

	# For outer tools
	for my $k (qw/exe_samtools exe_hisat2 pl_runHisat2_with2pass/) {
		defined $opts{$k} and $gg{$k} = $opts{$k}; 
	}

	for my $k (qw/pref wrk_dir/) {
		defined $opts{$k} and $gg{$k} = $opts{$k}; 
	}

	$gg{'inBam'} = "$gg{'wrk_dir'}/$gg{'pref'}_comb.nonTrans.bam"; 
	defined $opts{'inBam'} and $gg{'inBam'} = $opts{'inBam'}; 
	$gg{'outBam'} = "$gg{'wrk_dir'}/$gg{'pref'}_comb.nonTransRmdup.bam"; 
	defined $opts{'outBam'} and $gg{'outBam'} = $opts{'outBam'}; 

	return; 
}# applyOpt() 

sub setGlob {
	$gg{'pref'}    = 'pref'; 
	$gg{'wrk_dir'} = "./"; 
	$gg{'inBam'}   = "$gg{'wrk_dir'}/$gg{'pref'}_comb.nonTrans.bam"; 
	$gg{'outBam'}  = "$gg{'wrk_dir'}/$gg{'pref'}_comb.nonTransRmdup.bam"; 

	$gg{'exe_samtools'} = 'samtools'; 

$gg{'help_txt'} = <<"HH"; 
################################################################################
# perl $0   -pref $gg{'pref'} -wrk_dir $gg{'wrk_dir'}
#
#   -pref       [$gg{'pref'}] Output prefix 
#     -inBam    [$gg{'inBam'}]  Not necessary. 
#     -outBam   [$gg{'outBam'}] Not necessary. 
#   -wrk_dir    [$gg{'wrk_dir'}] Directory storing result files. 
#
#   -exe_samtools      [$gg{'exe_samtools'}]
################################################################################
HH
	return; 
}# setGlob() 


