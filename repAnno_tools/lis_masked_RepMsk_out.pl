#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and die "perl $0 query.fa.kl in_RepMsk.out\n"; 
my $klF = shift; 

my $min_ident = 80; 
my $min_cov = 90; 

open KL,'<',"$klF" or die;
my (@ids, %len); 
while (<KL>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[0] eq 'key' and next; 
	push(@ids, $ta[0]); 
	$len{$ta[0]} = $ta[1]; 
}
close KL; 

my %blks; 
while (<>) {
	unless ( m/^\s*\d+/ ) {
		# print; 
		next; 
	}
	chomp; 
	s/^\s+//; s/\s+$//; 
	my @ta = split(/\s+/, $_); 
	$ta[1] <= 100-$min_ident or next; 
	$ta[5] <= $ta[6] or die "$_\n"; 
	push(@{$blks{$ta[4]}}, [$ta[5], $ta[6]]); 
}

my %covL; 
for (keys %blks) {
	my @new_blk; 
	for my $tr1 ( sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$blks{$_}} ) {
		if ( scalar(@new_blk) > 0 ) {
			if ( $new_blk[-1][1]+1 >= $tr1->[0] ) {
				$new_blk[-1][1] < $tr1->[1] and $new_blk[-1][1] = $tr1->[1]; 
			} else {
				push(@new_blk, [@$tr1]); 
			}
		} else {
			push(@new_blk, [@$tr1]); 
		}
	}
	for my $tr2 ( @new_blk ) {
		$covL{$_} += ( $tr2->[1]-$tr2->[0]+1 ); 
	}
}

for (@ids) {
	( defined $covL{$_} and $covL{$_} >= $min_cov/100 * $len{$_} ) or next; 
	$covL{$_} > $len{$_} and die "Longer: $_\n"; 
	m/^(RR\d+)_/ or die "Name: $_\n"; 
	print "$_\t$1_\t$1\t$covL{$_}\t$len{$_}\n"; 
}
