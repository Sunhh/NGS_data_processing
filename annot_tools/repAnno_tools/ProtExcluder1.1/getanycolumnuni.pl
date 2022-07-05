#! /usr/bin/perl -w

$usage = "getanycolumnuni.pl file column wanted\n";

# to get a list from a given column in a text file
# redundancy is excluded if the same items are next to each other

if (@ARGV < 2) {die "$usage";}

open(MSP, "$ARGV[0]") || die "Can not open the input MSP file $ARGV[0]\n$usage";

$lquery = "";

while (<MSP>) {
    @line = split;
    $i = $ARGV[1] -1;
    if ($line[$i] ne $lquery) {
	printf "%s\n",$line[$i];
    }
    $lquery = $line[$i];
}
close MSP;

