#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"taxa_list:s", 
); 

my $help_txt = <<HH; 

perl $0 all_orthomcl.out -taxa_list taxa_list > cafe_input.tab 

Format of -tax_list : tax1 \\n tax2 \\n tax3 ...

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
defined $opts{'taxa_list'} or &LogInforSunhh::usage($help_txt); 

my %need_taxa = %{&load_taxa_lis( $opts{'taxa_list'} )}; 
my @need_taxa_arr = sort { $need_taxa{$a} <=> $need_taxa{$b} } keys %need_taxa; 
my $tax_num = scalar(keys %need_taxa); 

print STDOUT join("\t", qw/Description ID/, @need_taxa_arr)."\n"; 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my @tb = split(/\s+/, $ta[1]); 
	my %cnt; 
	for my $tc (@tb) {
		$tc =~ m/^\s*$/ and next; 
		$tc =~ m/^(\S+)\((\S+)\)$/ or &stopErr( "tc=[$tc]\n" ); 
		my ($gid, $taxID) = ($1, $2); 
		$cnt{$taxID}{$gid} ++; 
	}
	my $grpID = &grpID( $ta[0] ); 
	my $grpDesc = $grpID; 
	print STDOUT "$grpDesc\t$grpID"; 
	for my $tk ( @need_taxa_arr ) {
		my $num = ( defined $cnt{$tk} ) ? scalar(keys %{$cnt{$tk}}) : 0 ; 
		print STDOUT "\t$num"; 
	}
	print STDOUT "\n"; 
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

sub load_taxa_lis {
	my %lis; 
	open F,'<',"$_[0]" or die "$!\n"; 
	while (<F>) {
		chomp; 
		my @ta= split(/\t/, $_); 
		$lis{$ta[0]} //= $.; 
	}
	close F; 
	return(\%lis); 
}

