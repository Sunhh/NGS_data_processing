use strict; 
use warnings; 

-t and !@ARGV and die "perl $0 in_result_list\n"; 

-d "all_results" or mkdir("all_results/"); 
my %h; 
my $suff; 
while (<>) { 
	chomp; 
	m!^\./structure_(\d+)/structure_(K\d+)! or die "$_\n"; 
	$suff = "${1}_${2}"; 
	my $bn=$_; 
	$bn=~s!^.+/!!; 
	my $new_f = "${suff}_${bn}"; 
	defined $h{$new_f} and die "$new_f\n"; $h{$new_f} = 1; 
	system "cp -p $_ all_results/$new_f";  
} 

