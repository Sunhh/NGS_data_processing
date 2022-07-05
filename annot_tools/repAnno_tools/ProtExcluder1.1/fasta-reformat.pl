#! /usr/bin/perl -w

$usage = "fasta-reformat.pl input-fasta-file number-of-positions-per-line\n";

# to reformat fasta file so that each line containing number of letters given by the user

if (@ARGV < 2) {die $usage;}
if ($ARGV[1] < 1) {die $usage;}

open(FA, "$ARGV[0]") || die $usage;

$seq = "";
while (<FA>) {
    if (/>\s*(.+)/) {
	if ($seq) {
	    @sym = split(//, $seq);
	    $ct = 0;
	    foreach $sym (@sym) {
		print $sym;
		$ct ++;
		if ( !($ct%$ARGV[1]) ) {print "\n";}
	    }
	    if ($ct%$ARGV[1]) {print "\n";}
	}
	printf ">%s\n", $1;
	$seq = "";
    } else {
	chomp;
	$seq .= $_;
    }
}
close FA;

@sym = split(//, $seq);
$ct = 0;
foreach $sym (@sym) {
    print $sym;
    $ct ++;
    if ( !($ct%$ARGV[1]) ) {print "\n";}
}
if ($ct%$ARGV[1]) {print "\n";}
