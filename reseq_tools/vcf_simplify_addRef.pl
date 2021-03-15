#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and die "perl $0 refID sample.vcf\n"; 

my $id = shift; 

while (<>) {
	chomp; 
	m!^\s*##! and do { print "$_\n"; next; }; 
	my @ta = split(/\t/, $_); 
	if (m!^\s*#CHROM!) {
		print join("\t", @ta[0..8], $id, @ta[9..$#ta])."\n"; 
		next; 
	}
	$ta[7] = '.'; 
	my $ii; 
	if ($ta[8] =~ m!^GT(:|$)!) {
		$ii = 0; 
	} else {
		my @tc = split(/:/, $ta[8]); 
		for (my $i0=0; $i0<@tc; $i0++) {
			$tc[$i0] eq 'GT' and do { $ii = $i0; last; }; 
		}
	}
	defined $ii or do { &tsmsg("[Err][Wrn] bad FORMAT: [$ta[8]]: $_\n"); next; }; 
	$ta[8] = 'GT'; 
	for my $tb (@ta[9..$#ta]) {
		$tb = (split(/:/, $tb))[$ii]; 
	}
	print join("\t", @ta[0..8], '0/0', @ta[9..$#ta])."\n"; 
}

