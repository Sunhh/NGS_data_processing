#!/usr/bin/env perl 
use strict; 
use warnings; 
use mathSunhh; 
my $ms = mathSunhh->new(); 

-t and !@ARGV and die "perl $0 P1R02scaf.nt_bn6.MiORGN_join 1>sep_ex_lis 2>joined_ex_lis\n"; 

my %merged_blk; # Blocks for Ex class 
my %merged_in_blk; # Blocks for In class
my %merged_un_blk; # Blocks for Un-defined class that including In and Ex both. 
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

	my $is_in = -1; 
	my $is_ex = -1; 

	if ( $en > 0 ) {
		if ( $in == 0 ) {
			$is_ex = 1; 
		} elsif ( $in <= 3 ) {
			if ( $en >= $in+3 ) {
				$is_ex = 1; $is_in = 0; 
			} elsif ( $in >= $en+3 ) {
				$is_in = 1; $is_in = 0; 
			} else {
				$is_in = 0; $is_ex = 0; 
			}
		} else { 
			if ( $en >= $in*2 ) {
				$is_ex = 1; $is_in = 0; 
			} elsif ( $in >= $en*2 ) {
				$is_in = 1; $is_ex = 0; 
			} else {
				$is_in = 0; $is_ex = 0; 
			}
		}
	} elsif ( $in > 0 ) {
		$is_in = 1; $is_ex = 0; 
	} else {
		next; 
	}

	if ( $is_in == 1 ) {
		push(@{$merged_in_blk{$ta[0]}}, [$ta[2], $ta[3]]); 
	} elsif ( $is_ex == 1 ) {
		push(@{$merged_blk{$ta[0]}}, [$ta[2], $ta[3]]); 
		push(@{$lines{$ta[0]}}, [@ta]); 
		print "$_\n"; 
	} elsif ( $is_in == 0 and $is_ex == 0 ) {
		push(@{$merged_un_blk{$ta[0]}}, [$ta[2], $ta[3]]); 
	} elsif ( $is_in == -1 and $is_ex == -1 ) {
		next; 
	} else {
		die "$_\n"; 
	}

}

for my $tid (keys %merged_blk) {
	my $tar = $ms->mergeLocBlk( $merged_blk{$tid} ); 
	my $cnt_len = 0; 
	for my $tr (@$tar) {
		$cnt_len += ($tr->[1]-$tr->[0]+1); 
	}
	my $cnt_len_in = 0; 
	if (defined $merged_in_blk{$tid}) {
		my $tar = $ms->mergeLocBlk( $merged_in_blk{$tid} ); 
		for my $tr (@$tar) {
			$cnt_len_in += ($tr->[1]-$tr->[0]+1); 
		}
	}
	my $cnt_len_un = 0; 
	if (defined $merged_un_blk{$tid}) {
		my $tar = $ms->mergeLocBlk( $merged_un_blk{$tid} ); 
		for my $tr (@$tar) {
			$cnt_len_un += ($tr->[1]-$tr->[0]+1); 
		}
	}

	my @tb = @{$lines{$tid}[0]}; 
	print STDERR join("\t", $tid, $tb[1], $tar->[0][0], $tar->[-1][1], $cnt_len, @tb[5 .. $#tb], $cnt_len_in, $cnt_len_un)."\n"; 
}


