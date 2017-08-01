#!/usr/bin/perl -w
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"cpuN:i", 
	"skip_hlineN:i", # 1
	"sub_hlineN:i",  # 1
	"cmd_str:s", 
	"addOnce_outhl:i", 
	"help!", 
); 
$opts{'cpuN'} //= 20; 
$opts{'skip_hlineN'} //= 1; 
$opts{'sub_hlineN'}  //= 1; 
$opts{'cmd_str'} //= 'cat'; 
$opts{'addOnce_outhl'} //= $opts{'sub_hlineN'}; 

my $pm; 
$opts{'cpuN'} > 1 and $pm = &LogInforSunhh::get_pm( $opts{'cpuN'} ); 

my $help_txt = <<HH; 

perl $0 -cmd_str '$opts{'cmd_str'}' input_filename > out_file_name

-help

-cmd_str        [$opts{'cmd_str'}] This will be the command used to process input file. 
-cupN           [$opts{'cpuN'}]    Number of threads to be used. 

-sub_hlineN     [$opts{'sub_hlineN'}] Only when cpuN > 1; Number of lines in main file which will be put into each sub-file. 
-skip_hlineN    [$opts{'skip_hlineN'}] Only when cpuN > 1; Number of lines which will be skipped when merging sub-files. 
-addOnce_outhl  [value of Number_of_sub_hlineN] Only when cpuN > 1; Add once 'addOnce_outhl' lines to output once. 

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
-t and !@ARGV and &LogInforSunhh::usage($help_txt); 



my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
my $fh = \*STDIN; 
@ARGV and $fh = &openFH( $ARGV[0] , '<' ); 

if ($opts{'cpuN'} <= 1) {
	open O,'>',"$wrk_dir/base_0" or die; 
	while (<$fh>) {
		print O $_; 
	}
	close O; 
	@ARGV and close($fh); 
	my $cmd = $opts{'cmd_str'}; 
	if ( $cmd =~ s!__INPUT__!$wrk_dir/base_0!g ) {
	} else {
		$cmd = "$cmd $wrk_dir/base_0"; 
	}
	&exeCmd_1cmd("$cmd > $wrk_dir/base_0.o"); 
	open F,'<',"$wrk_dir/base_0.o" or die; 
	while (<F>) {
		print STDOUT $_; 
	}
	close F; 
	&fileSunhh::_rmtree($wrk_dir);
	exit(0); 
}

my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => $opts{'sub_hlineN'}, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" ); 
close($fh); 
for my $sfn (@sub_fn) {
	my $pid = $pm->start and next;
	my $cmd = $opts{'cmd_str'}; 
	if ( $cmd =~ s!__INPUT__!$sfn!g ) {
	} else {
		$cmd = "$cmd $sfn"; 
	}
	&exeCmd_1cmd("$cmd > $sfn.o"); 
	$pm->finish;
}
$pm->wait_all_children;

if ($opts{'addOnce_outhl'} > 0) {
	open F,'<',"$sub_fn[0].o" or die; 
	while (<F>) {
		$. > $opts{'addOnce_outhl'} and last; 
		print STDOUT $_; 
	}
	close F; 
}

for my $sfn (@sub_fn) {
	open F,'<',"$sfn.o" or die; 
	for (my $i=0; $i<$opts{'skip_hlineN'}; $i++) {
		<F>; 
	}
	while (<F>) {
		print STDOUT $_; 
	}
	close F; 
}
&fileSunhh::_rmtree($wrk_dir); 



