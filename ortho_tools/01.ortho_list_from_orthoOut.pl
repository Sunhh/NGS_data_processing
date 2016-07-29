#!/usr/bin/perl 
# 20160728 Add more description. 
#   Only keep Orthologous Groups (OGs) with 1-to-1 protein relationships. 
#   Only consider taxon defined in taxa_list. 
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"taxa_list:s", 
	"allow_miss!", 
); 

my $help_txt = <<HH; 
####################################################################################################
# perl $0 all_orthomcl.out > all_orthomcl.out.1to1_OGs
# 
#   Only OGs with at most one gene for each taxon will be kept. 
#   By default, I require all taxa existing in an OG. 
#
# -taxa_list   [in_taxa.list] Format: tax1.fa \\n tax2.fa \\n ...
# -allow_miss  [Boolean] Allow missing of taxa if given. 
# 
# -help
# 
####################################################################################################
HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

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
		defined $need_taxa{ $taxID } or next; 
		push(@gIDs, $gid); 
		defined $cnt{$taxID} and do { $is_bad = 1; last; }; 
		push(@{$cnt{$taxID}}, $gid); 
	}
	$is_bad == 1 and next; 
	$opts{'allow_miss'} or scalar(keys %cnt) == $tax_num or next; 
	print STDOUT "$ta[0]"; 
	for my $tk ( @need_taxa_arr ) {
		print STDOUT "\t" . join(" ;; ", @{$cnt{$tk}}); 
	}
	print STDOUT "\n"; 
}


sub load_taxa_lis {
	# Input format : taxID_1 \\n taxID_2 \\n ....
	my %lis; 
	open F,'<',"$_[0]" or die "$!\n"; 
	while (<F>) {
		chomp; 
		my @ta= split(/\t/, $_); 
		$ta[0] =~ s/\.fa$//; 
		$lis{$ta[0]} //= $.; 
	}
	close F; 
	return(\%lis); 
}

