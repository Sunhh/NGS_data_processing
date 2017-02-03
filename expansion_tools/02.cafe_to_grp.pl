#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"cafe_report:s", # Required as cafe output. 
	"pvalue_node:i", # Default -1 for all 
	"pvalue_thres:f", # Default 0.01 
); 

$opts{'pvalue_node'}  //= -1; 
$opts{'pvalue_thres'} //= 0.01; 
my $help_txt = <<HH; 

perl $0 -cafe_report input.cafe

-help               

-pvalue_node        [$opts{'pvalue_node'}] 
-pvalue_thres       [$opts{'pvalue_thres'}] 

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
defined $opts{'cafe_report'} or &LogInforSunhh::usage($help_txt); 

my $inFh = &openFH($opts{'cafe_report'}, '<'); 
my %info_hash; 
$info_hash{'is_table'} = 0; 
while (<$inFh>) {
	m/^\s*$/ and next; 
	s!\t$!!; 
	chomp; 
	my $out_pvalue = ''; 
	if ($info_hash{'is_table'} == 1) {
		my @ta = split(/\t/, $_); 
		if ( $opts{'pvalue_node'} < 0 ) {
			# Select overall p-value 
			$ta[2] <= $opts{'pvalue_thres'} and print STDOUT "$_\t$ta[2]\n"; 
			next; 
		}
		defined $info_hash{'slct_nodeIJ'} or &stopErr("[Err] No selected node information found.\n"); 
		# ((0.004210,0.000000),(0.010377,0.010010),(0.859910,0.560166)) 
		$ta[3] =~ s/^\(+(.+)\)+$/$1/ or &stopErr("[Err] Bad format of node_pvalue [$ta[3]]\n"); 
		my @p_pairs = split(/\)\s*,\s*\(/, $ta[3]); 
		for (my $i=0; $i<@p_pairs; $i++) {
			$i == $info_hash{'slct_nodeIJ'}[0] or next; 
			my @pp = split(/\s*,\s*/, $p_pairs[$i]); 
			$pp[0] =~ m/^\-$/ and next; 
			my $p_slct = $pp[ $info_hash{'slct_nodeIJ'}[1] ]; 
			$p_slct <= $opts{'pvalue_thres'} or next; 
			print STDOUT "$_\t$p_slct\n"; 
			last; 
		}
	} elsif ( m/^Tree:/ ) {
		next; 
	} elsif ( m/^Lambda:/ ) {
		next; 
	} elsif ( m/^Lambda tree:/ ) {
		next; 
	} elsif ( m/^# IDs of nodes:/ ) {
		next; 
	} elsif ( m/^# Output format for:\s*\'\s*Average Expansion/ ) {
		if (m/\'Branch\-specific P\-values\' = \(node ID, node ID\): \s*(\(.+\))\s*$/) {
			$opts{'pvalue_node'} < 0 and do { $info_hash{'slct_nodeIJ'} = [-1, -1]; next; }; 
			my $nodID_str = $1; 
			my @nodID_pairs = split(/\s+/, $nodID_str); 
			NODE_LIST: 
			for (my $i=0; $i<@nodID_pairs; $i++) {
				my $tb = $nodID_pairs[$i]; 
				$tb =~ m/^\((\d+),(\d+)\)$/ or &stopErr("[Err] Bad pair format: [$tb]\n"); 
				my @nn = ($1, $2); 
				for (my $j=0; $j<@nn; $j++) {
					if ($nn[$j] == $opts{'pvalue_node'}) {
						$info_hash{'slct_nodeIJ'} = [$i, $j]; 
						last NODE_LIST; 
					}
				}
			}
		} else {
			&stopErr("[Err] Failed to parse node ID: $_\n"); 
		}
		next; 
	} elsif ( m/^# Output format for\s+\'\s*Branch cutting/ ) {
		next; 
	} elsif ( m/^Average Expansion:/ ) {
		next; 
	} elsif ( m/^Expansion\s*:/ ) {
		next; 
	} elsif ( m/^Remain\s*:/ ) {
		next; 
	} elsif ( m/^Decrease:/ ) {
		next; 
	} elsif ( m/^\'ID\'/ ) {
		$info_hash{'is_table'} = 1; 
		s/\'//g; 
		print STDOUT "$_\tPvalue\n"; 
	} else {
		&stopErr("[Err] Some line I don't understand: $_\n"); 
	}
}
close($inFh); 


