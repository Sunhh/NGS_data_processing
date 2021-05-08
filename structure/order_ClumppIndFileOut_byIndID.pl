#!/usr/bin/perl
use strict; 
use warnings; 

!@ARGV and die "perl $0 order_list ClumppIndFile.output > ClumppIndFile.output.srt\n"; 

my $f_order = shift; 
my @id_ord; 
my %id; 
open F,'<',"$f_order" or die; 
while (<F>) {
	chomp; 
	my @ta= split(/\t/, $_); 
	push(@id_ord, $ta[0]); 
	defined $id{$ta[0]} and die "repeated ID [$ta[0]]\n"; 
	$id{$ta[0]} = $.; 
}
close F; 

my %lines; 
while (<>) {
	chomp; 
	my $raw_line = $_; 
	m!^\s*(\S+)! or die "$_\n"; 
	my $ii = $1; 
	defined $id{$ii} or do { warn "Skip missing ID [$ii]: $_\n"; next; }; 
	defined $lines{$ii} and die "repeat ID line [$_]\n"; 
	$lines{$ii} = $raw_line; 
}
for (@id_ord) {
	defined $lines{$_} or do { warn "Skip bad ID [$_]\n"; next; }; 
	print "$lines{$_}\n"; 
}

#   1        1   (0)      1 :  0.0000 1.0000
#   2        2   (2)      1 :  0.0000 1.0000
#   3        3   (0)      1 :  0.0000 1.0000
#   4        4   (0)      1 :  0.0000 1.0000
#   5        5   (0)      1 :  0.0000 1.0000
#   6        6   (0)      1 :  0.0000 1.0000


