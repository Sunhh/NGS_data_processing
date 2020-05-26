#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 

-t and !@ARGV and die "perl $0 R01/ct_R01.scafSeq.k27.counts > R01/ct_R01.scafSeq.k27.counts.tbl\n"; 

my (%seqs, $k, @klist); 
while (<>) {
	chomp; 
	if (m!^\>!) {
		m!^\>(\S+)! or die "$_\n"; 
		$k = $1; 
		defined $seqs{$k} and die "repeated k $k\n"; 
		push(@klist, $k); 
		next; 
	}
	s!^\s+!!; 
	s!\s+$!!; 
	$seqs{$k} .= " $_"; 
}
my $wdir = &fileSunhh::new_tmp_dir('create'=>1); 

my @ostat = qw/interval_mean interval_median interval_stdev MEAN MEDIAN MIN MAX NoNull/; 
print STDOUT join("\t", qw/SeqID/, @ostat)."\n"; 
for my $k1 (@klist) {
	$seqs{$k} =~ s!^\s+!!; 
	$seqs{$k} =~ s!\s+!\n!g; 
	&fileSunhh::write2file("$wdir/nn", "$seqs{$k}\n",'>'); 
	my %vstat = map { 
		split(/\t/, $_); 
	} `deal_table.pl $wdir/nn -col_stat 0 -col_stat_AsINS | deal_table.pl -transpose `; 
	print STDOUT join("\t", $k, @vstat{@ostat})."\n"; 
}

&fileSunhh::_rmtree($wdir); 

