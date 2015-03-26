#!/usr/bin/env perl 
use strict; 
use warnings; 
use mathSunhh; 
my $ms = mathSunhh->new(); 

-t and !@ARGV and die "perl $0 P1R02scaf.nt_bn6.MiORGN_join 1>sep_ex_lis 2>joined_ex_lis\n"; 

my %merged_blk; 
my %lines; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my ($in, $en) = (0,0); 
	if ($ta[6] =~ /In:([\d.]+)/) {
		$in = $1; 
	}
	if ($ta[6] =~ /Ex:([\d.]+)/) {
		$en = $1; 
	}

	$en > 0 or next; 
	
	if ($in == 0) {
		; 
	} elsif ($in <= 3) {
		$en >= $in+3 or next; 
	} else {
		$en >= $in*2 or next; 
	}
	push(@{$merged_blk{$ta[0]}}, [$ta[2], $ta[3]]); 
	push(@{$lines{$ta[0]}}, [@ta]); 
	print "$_\n"; 
}

for my $tid (keys %merged_blk) {
	my $tar = $ms->mergeLocBlk( $merged_blk{$tid} ); 
	my $cnt_len = 0; 
	for my $tr (@$tar) {
		$cnt_len += ($tr->[1]-$tr->[0]+1); 
	}
	my @tb = @{$lines{$tid}[0]}; 
	print STDERR join("\t", $tid, $tb[1], $tar->[0][0], $tar->[-1][1], $cnt_len, @tb[5 .. $#tb])."\n"; 
}


