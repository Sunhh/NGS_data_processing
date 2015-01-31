#! /usr/bin/perl

$usage = "rmlowcomfromMSP.pl mtca BFfile\n";

if (@ARGV < 2) {die $usage;}

# This is for finding the intems in file two  but not in file one

open(RM, "$ARGV[0]") || die "Cannot open $ARGV[0]";


while (<RM>) {
     if (/^>(\S+)\s+(\S+)\s*/) {
	 if ($1 ne $lpr && $2 ne $lseq) {
	 	$pr{$ct} = $1;
	$seq{$ct} = $2;
	$ct ++;
	 }
	 $lpr = $1;
	 $lseq = $2;
     }
}


open(MSP, "$ARGV[1]") || die "Cannot open $ARGV[1]";
while (<MSP>) {
         if (/^\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+(\S+)\s+\d+\s+\d+\s+\d+\s+(\S+)\s*/) {
	    	     if (&comparison) {
		 print;
	     }
	}
}
close MSP;

sub comparison  {
    foreach $key (keys %pr){
  if ($1 eq $seq{$key} && $2 eq $pr{$key} ) 	
  { return 1;}
      }
  }
