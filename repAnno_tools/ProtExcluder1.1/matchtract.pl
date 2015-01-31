#! /usr/bin/perl

$usage = "matchtract.pl blastx_output_file\n";

# to extract matched amino acids in blastx file

if (@ARGV < 1) {die "$usage";}
if (@ARGV > 1) {$score_cutoff = $ARGV[1];}
else {$score_cutoff = 0;}
if (@ARGV > 2) {$iden_cutoff = $ARGV[2];}
else  {$iden_cutoff = 0;}
open(BLT, "$ARGV[0]") || die "Can not open BLAST output $ARGV[0].\n$usage";

$score = -1;
while (<BLT>) {
    if (/^Query=\s+(\S+)/) {
	$query = $1;
    }
    elsif (/^>\s+(\S+)/) {
	$subject = $1;
	printf ">%s %s\n",$subject,$query;
    }
    elsif (/^Query\s+(\S+)/) {
	$take = 1;
    }
    elsif ($take) {
	print;
	$take = 0;
    }
}


