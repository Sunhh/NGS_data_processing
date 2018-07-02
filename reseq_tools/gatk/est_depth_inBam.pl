#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use mathSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"inBam:s", 
	"inBed:s", 
	"exe_samtools:s", 
); 

$opts{'exe_samtools'} //= 'samtools'; 

my $help_txt = <<HH; 
################################################################################
# perl $0 -inBam input_fix.bam -inBed input_cds.bed > input_fix.bam.depStat
#
# -exe_samtools     [samtools]
#
HH

defined $opts{'inBam'} or &LogInforSunhh::usage($help_txt); 


my $para_smtDep = ' -aa '; 
defined $opts{'inBed'} and $para_smtDep .= " -b $opts{'inBed'} "; 

open F,'-|',"samtools depth $para_smtDep $opts{'inBam'} | cut -f 3 " or die; 
my @Data; 
my $total = 0; 
while (<F>) {
	chomp; 
	push(@Data, $_); 
}
close F; 
$#Data == -1 and &stopErr("[Err] No input data found.\n"); 
&tsmsg("[Msg] Finish read in $opts{'inBam'}\n"); 
my @SortData = sort {$a<=>$b;} @Data; 
my ($min,$max,$mean) = ($SortData[0],$SortData[$#SortData],$total/($#SortData+1)); 
my $median;
if ($#SortData%2==0) {
	my $i = $#SortData/2; 
	$median = $SortData[$i]; 
} else {
	my $i1 = ($#SortData+1)/2; 
	my $i2 = $i1-1; 
	$median = ($SortData[$i1]+$SortData[$i2])/2; 
}
my $bh = &mathSunhh::ins_calc(\@Data, 0); 
my @use_tk   = qw/SUM MEAN MEDIAN min max COUNT interval_mean interval_median interval_stdev interval_low interval_high limit_low limit_high/; 
my @print_tk = qw/SUM MEAN MEDIAN MIN MAX Count interval_mean interval_median interval_stdev interval_low interval_high limit_low limit_high/; 
my @use_tv = @{$bh}{@use_tk}; 
print STDOUT join("\t", 'bamfile', @print_tk)."\n"; 
print STDOUT join("\t", $opts{'inBam'}, @use_tv)."\n"; 



