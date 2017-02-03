#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
); 

my $help_txt = <<HH; 

perl $0 all_orthomcl.out > all_orthomcl.out.tab

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

print STDOUT join("\t", qw/OrthoGrpID TaxID GeneID/)."\n"; 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my @tb = split(/\s+/, $ta[1]); 
	my $grpID = &grpID( $ta[0] ); 
	for my $tc (@tb) {
		$tc =~ m/^\s*$/ and next; 
		$tc =~ m/^(\S+)\((\S+)\)$/ or &stopErr( "tc=[$tc]\n" ); 
		my ($gid, $taxID) = ($1, $2); 
		print STDOUT join("\t", $grpID, $taxID, $gid)."\n"; 
	}
}

sub grpID {
	my $back = $_[0]; 
	if ($back =~ m/^(\S+)\s*\(\s*\d+\s+genes?\s*,\s*\d+\s*taxa\s*\)/) {
		$back = $1
	} else {
		$back =~ s/\s/_/g;
	}
	return $back; 
}
