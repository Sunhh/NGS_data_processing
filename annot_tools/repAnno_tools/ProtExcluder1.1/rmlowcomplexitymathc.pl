#! /usr/bin/perl

$usage = "gcalongenes.pl gc3 file aa number  minimum percent\n";


if (@ARGV < 3) {die $usage;}


open(MSP, "$ARGV[0]") || die $usage;
 
while (<MSP>){
    @line = split;
    $total = 0;
    $i = 3;
    $j = $i + $ARGV[1];
while ($i < $j) {
    $total = $total + $line[$i];
	$i ++;
    }
    if ($line[2]) {
    $rate = $total*100/$line[2];
    }
    if ($rate < $ARGV[2]) {
	print;
    }

}

close MSP;


