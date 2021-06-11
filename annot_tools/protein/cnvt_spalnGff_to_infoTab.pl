#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 

!@ARGV and die "perl $0 prot.key_len melon_v3.5raw.pep.fas.spaln.filter2.fmt.gff3 > melon_v3.5raw.pep.fas.spaln.filter2.fmt.gff3_info\n"; 

my $fnKL = shift; 

my $fh1 = &openFH($fnKL, '<'); 
my %k2l; 
while (<$fh1>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$k2l{$ta[0]} = $ta[1]; 
}
close($fh1); 

print STDOUT join("\t", qw/ScfID  CovByScf  CovS2  CovE2  CovPerc/)."\n"; 
my %h; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[2] eq 'match_part' or next; 
	$ta[8] =~ m!Parent=([^\s;]+);Target=([^\s;]+)\s+(\d+)\s+(\d+)\s+([+-])\s*(?:;|$)! or die "$_\n"; 
	my ($pid, $tid, $ts, $te, $tstr) = ($1, $2, $3, $4, $5); 
	defined $k2l{$tid} or die "Failed to know length of [$tid]\n"; 
	$h{$pid}{'refID'} //= $ta[0]; 
	$h{$pid}{'s'} //= $ta[3]; $h{$pid}{'s'} > $ta[3] and $h{$pid}{'s'} = $ta[3]; 
	$h{$pid}{'e'} //= $ta[4]; $h{$pid}{'e'} < $ta[4] and $h{$pid}{'e'} = $ta[4]; 
	$h{$pid}{'sum'} += ($te-$ts+1); 
	$h{$pid}{'order'} //= $.; 
	$h{$pid}{'tid'} //= $tid; 
}
for my $pid (sort { $h{$a}{'order'} <=> $h{$b}{'order'} } keys %h) {
	print STDOUT join("\t", $h{$pid}{'tid'}, $h{$pid}{'refID'}, $h{$pid}{'s'}, $h{$pid}{'e'}, sprintf("%.2f", 100*$h{$pid}{'sum'}/$k2l{$h{$pid}{'tid'}}))."\n"; 
}



