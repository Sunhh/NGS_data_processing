#!/usr/bin/perl
use strict; 
use warnings; 

!@ARGV and die "perl $0 in.win.fst\n"; 

my $f = shift; 

system "awk ' NR > 1 && ( \$4 >= 50 && \$6 >= 0 ) ' $f | deal_table.pl -col_stat 5 -col_stat_AsINS | deal_table.pl -transpose > $f.statMean"; 

