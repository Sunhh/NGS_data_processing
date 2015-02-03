#! /usr/bin/perl -w
 
$usage = "mergeunmatchedregion.pl seqfile\n";

# to merge multiple pieces from a single sequence into one piece
 
open(DB, "$ARGV[0]") || die "Can not open the seqfile $ARGV[0]\n$usage";
 
$lTE = "";
while (<DB>) {
    if (/^>(\S+)\D\d+-\d+\s*(.*)$/){
	if ($1 ne $lTE) {
	    printf ">%s\t %s\n",$1, $2;
	}
	$lTE = $1;
    } else {print;}
}

close DB;
