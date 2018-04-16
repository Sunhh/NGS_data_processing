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
	"noBase0!", 
	"dvdInput:i", 
	"help!", 
); 
$opts{'cpuN'} //= 20; 
$opts{'skip_hlineN'} //= 1; 
$opts{'sub_hlineN'}  //= 1; 
$opts{'cmd_str'} //= 'cat'; 
$opts{'addOnce_outhl'} //= $opts{'sub_hlineN'}; 
$opts{'dvdInput'} //= 1; 

my $pm; 
$opts{'cpuN'} > 1 and $pm = &LogInforSunhh::get_pm( $opts{'cpuN'} ); 

my $help_txt = <<HH; 

perl $0 -cmd_str '$opts{'cmd_str'}' input_filename > out_file_name

-help

-cmd_str        [$opts{'cmd_str'}] This will be the command used to process input file. 
-cpuN           [$opts{'cpuN'}]    Number of threads to be used. 

-sub_hlineN     [$opts{'sub_hlineN'}] Only when cpuN > 1; Number of lines in main file which will be put into each sub-file. 
-skip_hlineN    [$opts{'skip_hlineN'}] Only when cpuN > 1; Number of lines which will be skipped from sub-output-file when merging. 
-addOnce_outhl  [value of Number_of_sub_hlineN] Only when cpuN > 1; Add once 'addOnce_outhl' lines from first sub-output to merged_output before any sub-output-files merging. 

-noBase0        [Boolean] If given, won't use base_0 file as temporary file. 
-dvdInput       [$opts{'dvdInput'}] Divid the input into numbers of sub-files, and process each with multi-threading. 

"__INPUT__" could be in the cmd_str as a replace of input file in command. 

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
-t and !@ARGV and &LogInforSunhh::usage($help_txt); 



my $fh = \*STDIN; 
@ARGV and $fh = &openFH( $ARGV[0] , '<' ); 


if ($opts{'dvdInput'} > 1) {
	my $wrk_d0 = &fileSunhh::new_tmp_dir( 'create' => 1 ); 

	my @sub_f0; 
	if ($opts{'noBase0'}) {
		@sub_f0= &fileSunhh::dvd_file( $fh, $opts{'dvdInput'}, 'keep_order' => 1, 'with_header' => $opts{'sub_hlineN'}, 'sub_pref' => "$wrk_d0/sub_", 'tmpFile' => "" ); 
	} else {
		@sub_f0= &fileSunhh::dvd_file( $fh, $opts{'dvdInput'}, 'keep_order' => 1, 'with_header' => $opts{'sub_hlineN'}, 'sub_pref' => "$wrk_d0/sub_", 'tmpFile' => "$wrk_d0/base_0" ); 
	}
	close($fh); 

	my $has_addOnce = 0; 
	for my $f0 (@sub_f0) {
		my $fh0 = &openFH($f0, '<'); 
		my $fh1 = &openFH("$f0.o", '>'); 
		&proc_input($fh0, $fh1); 
		close($fh1); 
		if ($opts{'addOnce_outhl'} > 0 and $has_addOnce == 0) {
			open F,'<',"$f0.o" or die "$!\n"; 
			while (<F>) {
				$. > $opts{'addOnce_outhl'} and last; 
				print STDOUT $_; 
			}
			close F; 
			$has_addOnce = 1; 
		}
		open F,'<',"$f0.o" or die; 
		for (my $i=0; $i<$opts{'skip_hlineN'}; $i++) {
			<F>; 
		}
		while (<F>) {
			print STDOUT $_; 
		}
		close F; 
	}
	&fileSunhh::_rmtree($wrk_d0); 

} else {
	&proc_input($fh, \*STDOUT); 
}


##################################

sub proc_input {
	my ($i_fh, $o_fh) = @_; 
	my $i_wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 

	if ($opts{'cpuN'} <= 1) {
		open O,'>',"$i_wrk_dir/base_0" or die "$!\n"; 
		while (<$i_fh>) {
			print O $_; 
		}
		close O; 
		close($i_fh); 
		my $cmd = $opts{'cmd_str'}; 
		if ( $cmd =~ s!__INPUT__!$i_wrk_dir/base_0!g ) {
		} else {
			$cmd = "$cmd $i_wrk_dir/base_0"; 
		}
		&exeCmd_1cmd("$cmd > $i_wrk_dir/base_0.o"); 
		open F,'<',"$i_wrk_dir/base_0.o" or die; 
		while (<F>) {
			print {$o_fh} $_; 
		}
		close F; 
		&fileSunhh::_rmtree($i_wrk_dir); 
		return(); 
	}

	my @i_sub_fn; 
	if ($opts{'noBase0'}) {
		@i_sub_fn= &fileSunhh::dvd_file( $i_fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => $opts{'sub_hlineN'}, 'sub_pref' => "$i_wrk_dir/sub_", 'tmpFile' => "" ); 
	} else {
		@i_sub_fn= &fileSunhh::dvd_file( $i_fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => $opts{'sub_hlineN'}, 'sub_pref' => "$i_wrk_dir/sub_", 'tmpFile' => "$i_wrk_dir/base_0" ); 
	}
	close($i_fh); 

	for my $i_sfn (@i_sub_fn) {
		my $pid = $pm->start and next;
		my $cmd = $opts{'cmd_str'}; 
		if ( $cmd =~ s!__INPUT__!$i_sfn!g ) {
		} else {
			$cmd = "$cmd $i_sfn"; 
		}
		&exeCmd_1cmd("$cmd > $i_sfn.o"); 
		$pm->finish;
	}
	$pm->wait_all_children;

	if ($opts{'addOnce_outhl'} > 0) {
		open F,'<',"$i_sub_fn[0].o" or die; 
		while (<F>) {
			$. > $opts{'addOnce_outhl'} and last; 
			print {$o_fh} $_; 
		}
		close F; 
	}

	for my $i_sfn (@i_sub_fn) {
	open F,'<',"$i_sfn.o" or die; 
		for (my $i=0; $i<$opts{'skip_hlineN'}; $i++) {
			<F>; 
		}
		while (<F>) {
			print {$o_fh} $_; 
		}
		close F; 
	}
	&fileSunhh::_rmtree($i_wrk_dir); 
}# proc_input() 



