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


my %need_taxa = %{&load_taxa_lis( $opts{'taxa_list'} )}; 
my @need_taxa_arr = sort { $need_taxa{$a} <=> $need_taxa{$b} } keys %need_taxa; 
my $tax_num = scalar(keys %need_taxa); 

print STDOUT join("\t", qw/GrpInfor/, @need_taxa_arr)."\n"; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my @tb = split(/\s+/, $ta[1]); 
	my %cnt; 
	my $is_bad = 0; 
	my @gIDs; 
	for my $tc (@tb) {
		$tc =~ m/^\s*$/ and next; 
		$tc =~ m/^(\S+)\((\S+)\)$/ or die "tc=[$tc]\n"; 
		my ($gid, $taxID) = ($1, $2); 
		push(@gIDs, $gid); 
		defined $cnt{$taxID} and do { $is_bad = 1; last; }; 
		push(@{$cnt{$taxID}}, $gid); 
	}
	$is_bad == 1 and next; 
	scalar(keys %cnt) == $tax_num or next; 
	print STDOUT "$ta[0]"; 
	for my $tk ( @need_taxa_arr ) {
		print STDOUT "\t" . join(" ;; ", @{$cnt{$tk}}); 
	}
	print STDOUT "\n"; 
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

