#!/usr/bin/env perl
use strict; 
use warnings; 
use LogInforSunhh; 
use Parallel::ForkManager; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"cpuN:s", "beginLn:i", "endLn:i", "nprocF:s", 
	"grpLn:i", 
	"wait_sec:i", 
	"rmComment!", 
); 

sub usage {
	print STDOUT <<HH;
################################################################################
# perl $0 input_cmd_list
# 
# -beginLn       [0]
# -endLn         [0]
# -cpuN          [1]
# -nprocF        [Nproc]
#
# -grpLn         [1] Do be CAREFUL when using this!!! 
#
# -rmComment     [Boolean]
# -wait_sec      [0]
################################################################################
HH
	exit 1; 
}

$opts{'beginLn'} = $opts{'beginLn'} // 0; 
$opts{'endLn'} = $opts{'endLn'} // 0; 
$opts{'cpuN'} = $opts{'cpuN'} // 0; 
$opts{'nprocF'} = $opts{'nprocF'} // 'Nproc'; 
$opts{'grpLn'} //= 1; 
$opts{'grpLn'} >= 1 or &stopErr("[Err] -grpLn cannot be smaller than 1.\n"); 
$opts{'wait_sec'} //= 0; 

-t and !@ARGV and &usage(); 
defined $opts{'help'} and &usage(); 

my @joblist; 
while (<>) {
	chomp; 
	push(@joblist, $_); 
}
@joblist = &linearize(\@joblist); 
$opts{'rmComment'} and @joblist = &remove_comment(\@joblist); 
my $ttlN = scalar(@joblist); 
&tsmsg("[Rec] Total $ttlN lines.\n"); 
$opts{'beginLn'} <= 0 and $opts{'beginLn'} = 1; 
$opts{'endLn'} > $ttlN and $opts{'endLn'} = $ttlN; 
$opts{'endLn'} <= 0 and $opts{'endLn'} = $ttlN; 
&tsmsg("[Rec] Begin/End line modified to [$opts{'beginLn'} , $opts{'endLn'}]\n"); 

my %jobDone; 

my $MAX_PROCESSES = $opts{'cpuN'} ; # Sometimes $parm{'cpuN'} - 1 may be better.
my $pm = new Parallel::ForkManager($MAX_PROCESSES); 
for (my $i=$opts{'beginLn'}-1; $i<$opts{'endLn'}; $i+=$opts{'grpLn'}) {
	$MAX_PROCESSES = &change_procN($pm, $opts{'nprocF'}, $MAX_PROCESSES); 
	if ($opts{'wait_sec'} > 0) {
		sleep($opts{'wait_sec'}); 
	}
	my $pid = $pm->start and next; 
	my $i1 = $i+1; 
	for (my $j=0; $j<$opts{'grpLn'}; $j++) {
		my $k = $i+$j; 
		$k < $opts{'endLn'} or last; 
		my $cmdLn = $joblist[$k]; 
		if ( $cmdLn =~ m!^\s*$! ) {
			&tsmsg("[Wrn] Skip [$i1+$j] empty command : $cmdLn\n"); 
		} elsif ( $cmdLn =~ m!^\s*\#! ) {
			&tsmsg("[Wrn] Skip [$i1+$j] commented command : $cmdLn\n"); 
		} else {
			&tsmsg("[Msg] Running [$i1+$j] command : $cmdLn\n"); 
			&exeCmd_1cmd($cmdLn); 
		}
	}
	$pm->finish; 
	$MAX_PROCESSES = &change_procN($pm, $opts{'nprocF'}, $MAX_PROCESSES); 
}
$pm->wait_all_children; 

&tsmsg("[Rec] All commands over.\n"); 

################################################################################
# Sub-routines 
################################################################################


sub change_procN {
	my ($pm, $nprocF, $prev_maxP) = @_; 
	-e $nprocF or return $prev_maxP; 
	open F,'<',"$nprocF" or &stopErr("[Err] Failed to open [$nprocF].\n"); 
	my $new_maxP = <F>; 
	chomp($new_maxP); 
	$new_maxP = (split(/\s+/, $new_maxP))[0]; 
	$new_maxP = int($new_maxP); 
	close F; 
	if ($new_maxP > 0 and $new_maxP != $prev_maxP) {
		$pm->set_max_procs($new_maxP); 
		&tsmsg("[Rec] Changing MAX_PROCESSES from $prev_maxP to $new_maxP\n"); 
	}
	return $new_maxP; 
}# change_procN()


sub linearize {
	my ($in_aref) = @_; 
	my @back; 
	my $is_cont = 0; 
	for my $a1 (@$in_aref) {
		if ($is_cont == 1) {
			$back[-1] .= "$a1"; 
			$is_cont = 0; 
		} else {
			push(@back, $a1); 
		}
		$back[-1] !~ m!^\s*#! and $back[-1] =~ s!(\s)\\$!$1! and $is_cont = 1; 
	}
	if ( scalar(@back) < scalar(@$in_aref) ) {
		@back = &linearize(\@back); 
	}
	return (@back); 
}# linearize() 

sub remove_comment {
	my ($in_aref) = @_; 
	my $cnt = 0; 
	my @back; 
	for my $a1 (@$in_aref) {
		$a1 =~ m!^\s*(#|$)! and do { &tsmsg("[Wrn] Skip unusable line [$a1]\n"); next; }; 
		push(@back, $a1); 
	}
	return(@back); 
}# remove_comment() 


