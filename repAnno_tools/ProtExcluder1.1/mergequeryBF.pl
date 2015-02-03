#! /usr/bin/perl
 
$usage = "mergequeryBF.pl BF maximum gap to merge\n";

# to merge matched region in query if they are within given distance


if (@ARGV < 2) {die "$usage";}

`sort -k 6,6 -k 3,3n $ARGV[0] > $ARGV[0]s`;

open(MSP, "$ARGV[0]s") || die $usage;

while (<MSP>) {
  
if   (/^\s*\d+\s+\d+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+\d+\s+\d+\s+\d+\s+(\S+)\s*/)  {
    if ($4 ne $lTE) {
	printf "%-30s %06d %06d %06d %s\n", $lTE, $start,$end,$llen,$lsubj;
	$start = $1;
	$end = $2;
    }
    elsif (($1 - $end) > $ARGV[1]) {
	printf "%-30s %06d %06d %06d %s\n", $lTE, $start,$end,$llen,$lsubj;
	$start = $1;
	$end = $2;
    }
    elsif ($2 > $end) {
	$end = $2;
    }
    $llen = $3;
    $lTE = $4;
    $lsubj = $5;
    }
}

close MSP;

printf "%-30s %06d %06d %06d %s\n", $lTE, $start,$end,$llen,$lsubj;

