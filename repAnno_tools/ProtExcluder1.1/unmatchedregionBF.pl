#! /usr/bin/perl
 
$usage = "unmatchedregionBF.pl BFm50 bpflanking\n";

# to extract the unmatched portion of the sequence

if (@ARGV < 2) {die "$usage";}

open(BF, "$ARGV[0]") || die $usage;

while (<BF>) {
  
if   (/^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+\S+\s*/)  {
    $head = $lend + $ARGV[1];
    if ($1 ne $lTE) {
		if ($head < $llen) {
		printf "1000 85 1 200 query %06d %06d %s\n", $head, $llen,$lTE;
	}
	if ($2 > $ARGV[1]) {
	printf "1000 85 1 200 query 000001 %06d %s\n", ($2-$ARGV[1]),$1;
	}
    }
    elsif ($head < ($2 - $ARGV[1])) {
	printf "1000 85 1 200 query %06d %06d %s\n", $head, ($2-$ARGV[1]),$1;
    }
    $lTE = $1;
    $lend = $3;
    $llen = $4;
    }
}

close BF;

$head = $lend + $ARGV[1];
if ($head < $llen) {
    printf "1000 85 1 200 query %06d %06d %s\n", $head, $llen,$lTE;
}
