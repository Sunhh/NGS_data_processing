#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 

!@ARGV and die "perl $0 min_major_AN in.tab.lcn [cpuN]\n"; 

my $min_an = shift; 
my $in_fn = shift; 
my $cpuN = shift; 
$cpuN //= 30; 
my $pm = &LogInforSunhh::get_pm( $cpuN ); 

my $fh = &openFH( $in_fn, '<' ); 

my $header_txt = <$fh>; 
# print STDOUT $header_txt; 
print STDOUT join("\t", qw/chr pos/, $in_fn, "mj_AN", "any_AN")."\n"; 

my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
my @sub_fn = &fileSunhh::dvd_file( $fh, $cpuN, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" ); 
close($fh); 
for my $sfn (@sub_fn) {
	my $pid = $pm->start and next; 
	open F,'<',"$sfn" or &stopErr("[Err] Failed to open file [$sfn]\n"); 
	open O,'>',"$sfn.o" or &stopErr("[Err] Failed to open file [$sfn.o]\n"); 
	while (<F>) {
		chomp; 
		my @ta = &splitL("\t", $_); 
		$ta[6] =~ m!^([^\s_]+)_(\d+)! or die "line: $_\n"; 
		my ($mj_al, $mj_an) = ($1, $2); 
		my $any_an = 2 * ($ta[3]-$ta[4]); 
		if ( $any_an > 0 ) {
			$mj_an > 0.6 * $any_an or $mj_al = 'N'; 
		} else {
			$mj_al = 'N'; 
		}
		print O join("\t", @ta[0,1], $mj_al, $mj_an, $any_an)."\n"; 
	}
	close O; 
	close F; 
	$pm->finish; 
}
$pm->wait_all_children; 
for my $sfn (@sub_fn) {
	open F,'<',"$sfn.o" or die "Failed to open $sfn.o\n"; 
	while (<F>) {
		print STDOUT $_; 
	}
	close F; 
}
&fileSunhh::_rmtree($wrk_dir); 
