#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use mathSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"dist2join:i", # 1 
	"help!", 
); 
$opts{'dist2join'} //= 1; 

my $help_txt = <<HH; 

perl   $0    sep_loc    -dist2join $opts{'dist2join'}   > merged_loc 

-dist2join     [$opts{'dist2join'}]

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my $mat_obj = mathSunhh->new(); 

my %raw_blks; 
my %ord; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[1] =~ m/^[\d\.]+$/ or next; 
	$ord{$ta[0]} //= $.; 
	$ta[2] //= $ta[1]; 
	push(@{$raw_blks{$ta[0]}}, [$ta[1], $ta[2]]); 
}
my %merged_blks; 
print STDOUT join("\t", qw/chr start end/)."\n"; 
for my $tk (sort {$ord{$a} <=> $ord{$b}} keys %raw_blks) {
	$merged_blks{$tk} = $mat_obj->mergeLocBlk( $raw_blks{$tk}, 'dist2join'=>$opts{'dist2join'} ); 
	for my $tr1 (@{$merged_blks{$tk}}) {
		print STDOUT join("\t", $tk, $tr1->[0], $tr1->[1])."\n"; 
	}
}

