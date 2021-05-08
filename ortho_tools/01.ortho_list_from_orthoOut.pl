#!/usr/bin/perl 
# 20160728 Add more description. 
#   Only keep Orthologous Groups (OGs) with 1-to-1 protein relationships. 
#   Only consider taxon defined in taxa_list. 
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
use fileSunhh; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"taxa_list:s", 
	"allow_miss!", 
	"allow_mGene!", 
	"in_fmt:s", # orthomcl / orthofinder
); 

$opts{'in_fmt'} //= 'orthomcl'; 

my $help_txt = <<HH; 
####################################################################################################
# perl $0 all_orthomcl.out > all_orthomcl.out.1to1_OGs
# 
#   Only OGs with at most one gene for each taxon will be kept. 
#   By default, I require all taxa existing in an OG. 
#
# -taxa_list   [in_taxa.list] Format: tax1.fa \\n tax2.fa \\n ...
# -allow_miss  [Boolean] Allow missing of taxa if given. 
# -allow_mGene [Boolean] Allow multiple genes within each OG. 
#
# -in_fmt      [$opts{'in_fmt'}] orthomcl / orthofinder
# 
# -help
# 
####################################################################################################
HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

my %need_taxa = %{&load_taxa_lis( $opts{'taxa_list'} )}; 
my @need_taxa_arr = sort { $need_taxa{$a} <=> $need_taxa{$b} } keys %need_taxa; 
my $tax_num = scalar(keys %need_taxa); 

print STDOUT join("\t", qw/GrpInfor GrpID/, @need_taxa_arr)."\n"; 
if ( $opts{'in_fmt'} =~ m!^orthomcl$!i ) {
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
			unless ( $opts{'allow_mGene'} ) {
				defined $cnt{$taxID} and do { $is_bad = 1; last; }; 
			}
			push(@{$cnt{$taxID}}, $gid); 
		}
		$is_bad == 1 and next; 
		$opts{'allow_miss'} or scalar(keys %cnt) == $tax_num or next; 
		my $grpID = $ta[0]; 
		$grpID =~ s!\s*\(.*\)\s*(:\s*)?$!!; 
		print STDOUT "$ta[0]\t$grpID"; 
		for my $tk ( @need_taxa_arr ) {
			$cnt{$tk} //= ['NA']; 
			print STDOUT "\t" . join(" ;; ", @{$cnt{$tk}}); 
		}
		print STDOUT "\n"; 
	}

} elsif ( $opts{'in_fmt'} =~ m!^orthofinder$!i ) {
	my $h_l = <>; chomp($h_l); 
	my @h_a = &splitL("\t", $h_l); 
	my @h_i;
	for (my $i=0; $i<@need_taxa_arr; $i++) {
		my $tk; 
		for (my $j=1; $j<@h_a; $j++) {
			if ( $need_taxa_arr[$i] eq $h_a[$j] ) {
				$tk = $j; 
				last; 
			}
		}
		if ( defined $tk ){
			push(@h_i, $tk); 
		} else {
			&stopErr("[Err] No taxa [$need_taxa_arr[$i]] found in input.\n"); 
		}
	}
	# print STDOUT join("\t", $h_a[0], $h_a[0], @h_a[@h_i])."\n"; 

	while (<>) {
		chomp; s![^\S\t]+$!!; 
		my @ta = &splitL("\t", $_); 
		my %cnt; 
		my $is_bad = 0; 
		my @gIDs; 
		for my $j (@h_i) {
			unless ( $opts{'allow_miss'} ) {
				$ta[$j] =~ m!^\s*$! and do { $is_bad = 1; last; }; 
			}
			unless ( $opts{'allow_mGene'} ) {
				$ta[$j] =~ m!,! and do { $is_bad = 1; last; }; 
			}
			my @gid = grep { $_ !~ m/^\s*$/ } split(/,/, $ta[$j]) ; 
			if ( scalar(@gid) == 0 or (scalar(@gid) == 1 and $gid[0] =~ m!^\s*$!) ) {
				$gid[0] = ('NA'); 
			}
			push(@gIDs, join(" ;; ", @gid)); 
		}
		$is_bad == 1 and next; 
		print join("\t", $ta[0], $ta[0], @gIDs)."\n"; 
	}
} else {
	&stopErr("[Err] unknown -in_fmt [$opts{'in_fmt'}]\n"); 
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

