#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"in_orthomcl:s", 
	"in_name_list:s", 
); 

my $help_txt = <<HH; 

perl $0 -in_name_list xx.all_name_list   -in_orthomcl all_orthomcl.out > all_orthomcl.out.rawID 

HH

( defined $opts{'in_orthomcl'} and defined $opts{'in_name_list'} ) or &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my %id_orth2raw = %{ &load_name_list($opts{'in_name_list'}) }; 
my $fh_orth = &openFH($opts{'in_orthomcl'}, '<'); 
while (<$fh_orth>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my @tb = split(/\s+/, $ta[1]); 
	for my $tc (@tb) {
		$tc =~ m/^\s*$/ and next; 
		$tc =~ m/^(\S+)\((\S+)\)$/ or &stopErr( "tc=[$tc]\n" ); 
		my ($gid, $taxID) = ($1, $2); 
		my $new_gid = (defined $id_orth2raw{ $gid }) ? $id_orth2raw{ $gid } : $gid ; 
		$tc = "${new_gid}($taxID)"; 
	}
	print STDOUT join("\t", $ta[0], join(" ", @tb))."\n"; 
}
close ($fh_orth); 


# Not used here. 
sub grpID {
	my $back = $_[0]; 
	if ($back =~ m/^(\S+)\s*\(\s*\d+\s+genes?\s*,\s*\d+\s*taxa\s*\)/) {
		$back = $1
	} else {
		$back =~ s/\s/_/g;
	}
	return $back; 
}

sub load_name_list {
	# $_[0] : file name of in_name_list 
	my $fh = &openFH($_[0], '<'); 
	my %back; 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		$back{$ta[0]} = $ta[1]; 
	}
	close($fh); 
	return(\%back); 
}

