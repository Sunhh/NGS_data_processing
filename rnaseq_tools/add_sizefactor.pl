#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 

!@ARGV and die "perl $0 in.total.sizefactor in.cnt\n"; 

open F1,'<',"$ARGV[0]" or die; 
my %h; 
while (&wantLineC(\*F1)) {
	my @ta=&splitL("\t", $_); 
	$ta[0] eq 'sampleNames' and next; 
	# $ta[0] =~ s/^\S+_([FP][13]_[^_]+_rep\d+)$/$1/ or die "$_\n"; 
	$ta[0] =~ s/^\S+_([SM](?:FL|FR|LV|RT|SD|ST)(?:F1|P1|P3)_rep\d+)/$1/ or die "$ta[0]\n"; 
	$h{$1} = $ta[1]; 
}
close F1; 

open F2,'<',"$ARGV[1]" or die; 
while (<F2>) {
	print "$_"; 
	if ($. == 1) {
		chomp; 
		my @ta = split(/\t/, $_); 
		my @tb = ('sizeFactor'); 
		for (my $i=1; $i<@ta; $i++) {
			# $ta[$i] =~ s/^\S+_([FP][13]_[^_]+_rep\d+)$/$1/ or die "$ta[$i]\n"; 
			$ta[$i] =~ s/^\S+_([SM](?:FL|FR|LV|RT|SD|ST)(?:F1|P1|P3)_rep\d+)$/$1/ or die "$ta[$i]\n"; 
			defined $h{$ta[$i]} or die "$ta[$i]\n"; 
			$tb[$i] = $h{$ta[$i]}; 
		}
		print join("\t", @tb)."\n"; 
	}
}
close F2; 

