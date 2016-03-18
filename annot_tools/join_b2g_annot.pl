#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"gene_list:s", 
); 

my %hash; 
if (defined $opts{'gene_list'}) {
	my $fh = &openFH($opts{'gene_list'}, '<'); 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		defined $hash{$ta[0]} and next; 
		$hash{$ta[0]}[0] = $.; 
	}
	close($fh); 
}

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	defined $opts{'gene_list'} and !(defined $hash{$ta[0]}) and next; 
	$hash{$ta[0]}[0] //= $. ; 
	if ( $ta[1] =~ m/^GO:/ ) {
		push(@{$hash{$ta[0]}[1]}, $ta[1]); # GO ID 
	} elsif ( $ta[1] =~ m/^EC:/ ) {
		push(@{$hash{$ta[0]}[2]}, $ta[1]); # EC NO. 
	} else {
		die "$_\n"; 
	}
	if (defined $ta[2] and $ta[2] ne '') {
		push(@{$hash{$ta[0]}[3]}, $ta[2]); # B2G desc. 
	}
}

print join("\t", qw/GeneID GO GO_Num EC EC_Num Description Desc_Num/)."\n"; 
for my $tk (sort { $hash{$a}[0] <=> $hash{$b}[0] } keys %hash) {
	$hash{$tk}[1] //= []; 
	$hash{$tk}[2] //= []; 
	$hash{$tk}[3] //= []; 
	my $rmR_go = &rm_redund( $hash{$tk}[1] ); 
	my $rmR_ec = &rm_redund( $hash{$tk}[2] ); 
	my $rmR_desc = &rm_redund( $hash{$tk}[3] ); 
	my $n_go = scalar(@$rmR_go); 
	my $n_ec = scalar(@$rmR_ec); 
	my $n_desc = scalar(@$rmR_desc); 
	print join("\t", $tk, join(' ;; ', @$rmR_go), $n_go, join(' ;; ', @$rmR_ec), $n_ec, join(' ;; ', @$rmR_desc), $n_desc)."\n"; 
}

sub rm_redund {
	my @back; 
	my %h; 
	for (@{$_[0]}) {
		defined $h{$_} and next; 
		$h{$_} = 1; 
		push(@back, $_); 
	}
	return \@back; 
}


