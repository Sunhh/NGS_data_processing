#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"taxa_list:s", 
	"in_fmt:s", # orthomcl / orthofinder
); 

$opts{'in_fmt'} //= 'orthomcl'; 

my $help_txt = <<HH; 

perl $0 all_orthomcl.out -taxa_list taxa_list > cafe_input.tab 

-in_fmt     [$opts{'in_fmt'}] orthomcl / orthofinder

Format of -tax_list : tax1 \\n tax2 \\n tax3 ...

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
defined $opts{'taxa_list'} or &LogInforSunhh::usage($help_txt); 

my %need_taxa = %{&load_taxa_lis( $opts{'taxa_list'} )}; 
my @need_taxa_arr = sort { $need_taxa{$a} <=> $need_taxa{$b} } keys %need_taxa; 
my $tax_num = scalar(keys %need_taxa); 

print STDOUT join("\t", qw/Description ID/, @need_taxa_arr)."\n"; 

if ( $opts{'in_fmt'} =~ m/^orthomcl$/i ) {
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
} elsif ( $opts{'in_fmt'} =~ m/^orthofinder$/i ) {
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
		if (defined $tk) {
			push(@h_i, $tk); 
		} else {
			&stopErr("[Err] No taxa [$need_taxa_arr[$i]] found in input.\n"); 
		}
	}
	while (<>) {
		chomp; 
		my @ta = &splitL("\t", $_); 
		my $grpID = $ta[0]; 
		my $grpDesc = $grpID; 
		my @tb = @ta[@h_i]; 
		my @tc; 
		for (my $i=0; $i<@tb; $i++) {
			$tc[$i] = scalar( grep { $_ !~ m/^\s*$/ } split(/,/, $tb[$i]) ); 
		}
		print STDOUT join("\t", $grpID, $grpDesc, @tc)."\n"; 
	}
} else {
	&stopErr("[Err] Unknown format [$opts{'in_fmt'}]\n"); 
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

