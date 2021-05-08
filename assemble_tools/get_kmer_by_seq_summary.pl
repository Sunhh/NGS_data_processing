#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 

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
for my $v1 (@klist) {
	if (!defined $seqs{$v1}) {
		print STDOUT join("\t", $v1, (('0') x scalar(@ostat)))."\n"; 
		next; 
	}
	# &tsmsg("[Msg] Calculating for [$v1]\n"); 
	$seqs{$v1} =~ s!^\s+!!; 
	$seqs{$v1} =~ s!\s+!\n!g; 
	&fileSunhh::write2file("$wdir/nn", "$seqs{$v1}\n",'>'); 
	my %vstat = map { 
		chomp($_); 
		split(/\t/, $_); 
	} `deal_table.pl $wdir/nn -col_stat 0 -col_stat_AsINS | deal_table.pl -transpose `; 
	print STDOUT join("\t", $v1, @vstat{@ostat})."\n"; 
}

&fileSunhh::_rmtree($wdir); 

