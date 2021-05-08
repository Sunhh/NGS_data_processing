#!/usr/bin/perl
use strict; 
use warnings; 

print STDOUT <<TTT0; 

{
    "TaxaList":[
TTT0
my $cnt = 0; 
my @taxaArr; 
while (<>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	$cnt ++; 
	my $taxaTxt = <<TTT1; 
        {
            "name":"$ta[0]"
        }
TTT1
	chomp($taxaTxt); 
	push(@taxaArr, $taxaTxt); 
}
print STDOUT join(",\n", @taxaArr)."\n"; 
my $tailTxt = <<TTT2; 
    ]
}
TTT2
chomp($tailTxt); 
print STDOUT $tailTxt; 

