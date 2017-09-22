#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

# Depth file: 
# scaffold11      15      1
# scaffold11      16      1
# scaffold11      17      1
# scaffold11      18      1
# scaffold11      19      1

my $min_depth = 1; 
my $max_depth = -1; 

-t and !@ARGV and die "perl $0 in.bam.depth [minDepth maxDepth] > in.bam.depth.loc\nOr\ncat in.bam.depth | perl $0 [minDepth maxDepth] > in.bam.depth.loc\n"; 

if ( -t ) {
	scalar(@ARGV) >= 3 and ($max_depth) = splice( @ARGV, 2, 1 ); 
	scalar(@ARGV) >= 2 and ($min_depth) = splice( @ARGV, 1, 1 ); 
} else {
	scalar(@ARGV) >= 2 and ($max_depth) = splice( @ARGV, 1, 1 ); 
	scalar(@ARGV) >= 1 and ($min_depth) = splice( @ARGV, 0, 1 ); 
}

my @loci; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$min_depth > 0 and $ta[2] < $min_depth and next; 
	$max_depth > 0 and $ta[2] > $max_depth and next; 
	if (defined $loci[0]) {
		if ($loci[-1][0] eq $ta[0]) {
			if ( $loci[-1][2]+1 >= $ta[1] ) {
				$loci[-1][2] < $ta[1] and $loci[-1][2] = $ta[1]; 
			} else {
				push(@loci, [$ta[0], $ta[1], $ta[1]]); 
			}
		} else {
			push(@loci, [$ta[0], $ta[1], $ta[1]]); 
		}
	} else {
		push(@loci, [$ta[0], $ta[1], $ta[1]]); 
	}
}

for (@loci) {
	print STDOUT join("\t", @$_)."\n"; 
}

