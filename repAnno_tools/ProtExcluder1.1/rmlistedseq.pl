#! /usr/bin/perl -w

$usage = "rmlistedseq.pl namelist fasta file \n";

# This is for removing a subset of sequences from a seqfile

if (@ARGV < 2) {die $usage;}

open(RM, "$ARGV[0]") || die $usage;

$ct = 0;

while (<RM>) {
    if (/^>*(\S+)\s*/) {
        $seq{$ct} = $1;
	$ct ++;
    }
}
close RM;

open(FA, "$ARGV[1]") || die $usage;

while (<FA>) {
    if (/^>(\S+)\s*/) {
	if (&comparison) { 
	    $take = 0;
	}
	else {
	    $take = 1;
	}
    }
    if ($take) {
	print;
    }
}

close FA;




 
sub comparison {
 foreach $key (keys %seq){
     if ($1 eq $seq{$key}){ 
	 return 1;
     }
 }

}



