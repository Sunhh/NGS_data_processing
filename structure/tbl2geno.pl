use strict; 
use warnings; 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$. == 1 and $ta[1] !~ m/^\d+$/ and next; 
	$ta[1] =~ m/^\d+$/ or die "bad line: $_\n"; 
	print join("\t", @ta[0,1], join(" ", @ta[2..$#ta]))."\n"; 
}
