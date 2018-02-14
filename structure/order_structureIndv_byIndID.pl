#!/usr/bin/perl
use strict; 
use warnings; 

!@ARGV and die "perl $0 order_list P2_var_rep_3_f > P2_var_rep_3_f.srt\n"; 

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
my $is_Qh1 = 0; 
my $is_Qh2 = 0; 
my $done_Q = 0; 
my ($o_sec1, $o_sec2, $o_sec3); 
my %o_sec2_line; 
my $qNum = 0; 
while (<>) {
	chomp; 
	my $raw_line = $_; 
	if ( m!^Inferred ancestry of individuals! ) {
		$is_Qh1 = 1; 
		$o_sec1 .= "$_\n"; 
		next; 
	}
	if ( $is_Qh1 == 1 and m!^\s*Label\s+\(\%Miss\)\s+:\s+Inferred\s+clusters! ) {
		$is_Qh2 = 1; 
		$o_sec1 .= "$_\n"; 
		next; 
	}
	if ( $is_Qh2 == 1 and (m!(^\s*$)|^Estimated Allele Frequencies! or !(m!\s*(\d+)\s+\S+!))) {
		$is_Qh1 = $is_Qh2 = 0; 
		$done_Q = 1; 
		$o_sec3 .= "$_\n"; 
		next; 
	}
	if ($done_Q == 1) {
		$o_sec3 .= "$_\n"; 
	} elsif ( $is_Qh1 == 0 and $is_Qh2 == 0 ) {
		$o_sec1 .= "$_\n"; 
	} else {
		m!^(\s*)(\d+)(\s+\S.+:.+)$! or die "bad line:$_\n"; 
		my ($p1, $p2, $p3) = ($1, $2, $3); 
		defined $id{$p2} or die "Missing order for ID [$p2]:$_\n"; 
		defined $o_sec2_line{$p2} and die "repeat ID line [$_]\n"; 
		$qNum ++; 
		my $p2n = $id{$p2}; 
		$o_sec2_line{$p2} = [ $id{$p2}, "$p1$p2n$p3"]; 
	}
}
print $o_sec1; 
for my $k (sort {$o_sec2_line{$a}[0] <=> $o_sec2_line{$b}[0]} keys %o_sec2_line) {
	print "$o_sec2_line{$k}[1]\n"; 
}
print $o_sec3; 


#   1        1   (0)      1 :  0.0000 1.0000
#   2        2   (2)      1 :  0.0000 1.0000
#   3        3   (0)      1 :  0.0000 1.0000
#   4        4   (0)      1 :  0.0000 1.0000
#   5        5   (0)      1 :  0.0000 1.0000
#   6        6   (0)      1 :  0.0000 1.0000


