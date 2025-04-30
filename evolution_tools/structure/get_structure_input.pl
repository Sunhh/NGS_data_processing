#/usr/bin/perl
use strict;
use warnings;
## This program is to transform individual genotype data to a standard input file for the software "structure".
##Input file for this program is the genotype file for 50 Tibetans(all 22 auto chromosomes), if "genotype data with chimp" is used, it's OK, sequences beyond 50 are disgarded. 
# my %iupac = ("M" => "AC", 
#             "K" => "GT",   
#			 "Y" => "CT",
#			 "R" => "AG", 
#			 "W" => "AT", 
#			 "S" => "CG",
#			 "A" => "AA",
#			 "T" => "TT",
#			 "C" => "CC",
#			 "G" => "GG",
#			 "-" => "--",
#);
## Transform ATGC to 1234 and "-" to "-9" respectively. A->1 , T->2, C->3, G->4, - -> -9;
die "Usage:<input><list><output>\n" unless (@ARGV==3);

my %iupac = (
             "M" => [1,3],
			 "K" => [4,2],
			 "Y" => [3,2],
			 "R" => [1,4],
			 "W" => [1,2],
			 "S" => [3,4],
			 "A" => [1,1],
			 "T" => [2,2],
			 "C" => [3,3],
			 "G" => [4,4],
			 "-" => [-9,-9],
			 "N" => [-9, -9], 
);

open (IN, "<$ARGV[0]" )  || die "Can't open input file:[$ARGV[0]] $!\n";
open (INN, "<$ARGV[1]")   || die "Can't open input file:[$ARGV[1]] $!\n";
open (OUT, ">$ARGV[2]") || die "Can't open output file:[$ARGV[2]] $!\n";
#Store each line in an array, the first line is the distance between the two loci;
my @name_list =();
my $cout_num=0;
while(<INN>)
{
	chomp;
	my @inf=split ;
	$name_list[$cout_num]=$inf[-3];
	$cout_num++;
}
close INN;
print "sample num is $cout_num\n";

my @all_rows;
$all_rows[0] = "  -1";
for (my $i=1; $i<=$cout_num ; $i++) {
	$all_rows[$i] = ["$name_list[$i-1]", "$name_list[$i-1]",]; 
}

my $pre_chr = '';
my $pre_position;
while (<IN>) {
	chomp;
	next unless $_; #if an empty line is input, skip it;
	my @temp_line = split;
	$pre_chr  = $temp_line[0] unless $pre_chr; #For the first line, give $pre_chr a value;
	
	if ($temp_line[0] ne $pre_chr) {
		$all_rows[0] .= " -1";
		$pre_chr      = $temp_line[0]; #store the previous chr name;
		$pre_position = $temp_line[1]; # store the previous position;
	}elsif($pre_position){ #From the second time we read in a line, it goes this way.
        my $distance = $temp_line[1] - $pre_position; 
		$all_rows[0] .= " $distance";
	}
	$pre_position = $temp_line[1]; # store the position in order to calculate distrance between two adjacent loci;

	for (my $i = 2; $i <= $cout_num+1; $i++) {
defined $iupac{$temp_line[$i]} or die "|$temp_line[$i]|\n"; 
		$all_rows[$i-1]->[0] .= " $iupac{$temp_line[$i]}->[0]"; 
		$all_rows[$i-1]->[1] .= " $iupac{$temp_line[$i]}->[1]";
	}
}

my $count = @all_rows;
print OUT "$all_rows[0]\n";
for (my $j=1; $j<$count; $j++) {
	print OUT "$all_rows[$j]->[0]\n";
	print OUT "$all_rows[$j]->[1]\n";
}
close IN;
close OUT;
