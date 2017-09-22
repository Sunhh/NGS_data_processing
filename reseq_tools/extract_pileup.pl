#!/usr/bin/perl -w 
use strict; 

!@ARGV and die "perl $0 fsadf\n"; 

for (`ls *.pileup`) { 
	chomp; 
	print "Start $_;" . scalar(localtime()) . " \n"; 
	system "uniqComb.pl $_ -index ../basic_SNP_WM97toPI -col 0,1 -newCol 0,1 -exist > $_.basic"; 
	print "End $_;" . scalar(localtime()) . "\n"; 
} 
print "All extractions are over.\n";

