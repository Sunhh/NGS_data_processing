#!/usr/bin/perl
# 20150209 : Try to remove gene models who are overlapped by repeats. 
use strict; 
use warnings; 
use fileSunhh; 
use mathSunhh; 

use LogInforSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"debug!", 
	"idxGff:s", 
	"srcGff:s", 
	"ovl_ratio:f", 
); 

sub usage {
	print STDOUT <<HH; 
################################################################################
# perl $0 -srcGff gene_models.gff -idxGff repeat_loc.gff 
# 
# -ovl_ratio      [0.5] Maximum overlap ratio to keep a gene model. 
# 
# -help 
################################################################################
HH
	exit 1; 
}

defined $opts{'srcGff'} or &usage(); 
defined $opts{'idxGff'} or &usage(); 
$opts{'ovl_ratio'} = $opts{'ovl_ratio'} // 0.5; 
-t and !@ARGV and &usage(); 

my $srcFh = &openFH($opts{'srcGff'}, '<'); 
my $idxFh = &openFH($opts{'idxGff'}, '<'); 

# ##gff-version 3
# S401083_pilon   AUGUSTUS        match   113     6130    0.39    -       .       ID=0:g1.t1;
# S401083_pilon   AUGUSTUS        match_part      113     472     1       -       0       ID=0:g1.t1.c1;Parent=0:g1.t1;
# S401083_pilon   AUGUSTUS        match_part      595     807     0.99    -       0       ID=0:g1.t1.c2;Parent=0:g1.t1;
# S401083_pilon   AUGUSTUS        match_part      1269    1412    0.91    -       0       ID=0:g1.t1.c3;Parent=0:g1.t1;

my %rep_loc; 

while (<$idxFh>) {
	m/^#/ and next; 
	m/^\s*$/ and next; 
	m/^>/ and last; 
	chomp; 
	my @ta = split(/\t/, $_); 
	push(@{$rep_loc{$ta[0]}}, [@ta[3,4]]); 
}

%rep_loc = map { $_ => &merged_loc( $rep_loc{$_} );   } keys %rep_loc; 
sub merged_loc {
	my $ar = shift; 
	my @aa = sort { $a->[0]<=>$b->[0] || $a->[1] <=> $b->[1] } @$ar; 
	my @back; 
	for my $r1 (@aa) {
		my ($cs, $ce) = @$r1; 
		if ( scalar(@back) == 0 ) {
			@back = ([$cs,$ce]); 
		} else {
			if ( $back[-1][1] >= $cs-1 ) {
				$ce > $back[-1][1] and $back[-1][1] = $ce; 
			} else {
				push(@back, [$cs, $ce]); 
			}
		}
	}
	return \@back; 
}



my @src_geneL; 
while (<$srcFh>) {
	m/^#/ and next; 
	m/^\s*$/ and next; 
	m/^>/ and last; 
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($ta[2] eq 'mRNA') {
		if (scalar(@src_geneL) > 0) {
			&is_notOvl(\@src_geneL, \%rep_loc, $opts{'ovl_ratio'}) or &out_gL(\@src_geneL); 
		}
		@src_geneL = (); 
		push(@src_geneL, [[@ta], 'mRNA']); 
	} elsif ( $ta[2] eq 'CDS' ) {
		push(@src_geneL, [[@ta], 'CDS']); 
	} elsif ( $ta[2] =~ m/^(mRNA|exon|five_prime_UTR|three_prime_UTR)$/ ) {
		push(@src_geneL, [[@ta], $ta[2]]); 
	} else {
		&stopErr("[Err] Bad line: $_\n"); 
	}
}

if (scalar(@src_geneL) > 0) {
	&is_notOvl(\@src_geneL, \%rep_loc, $opts{'ovl_ratio'}) or &out_gL(\@src_geneL); 
	@src_geneL = (); 
}


sub is_notOvl {
	my ($ar, $hr, $ovl_ratio) = @_; 
	defined $ovl_ratio or &stopErr("[Err] Overlap ratio needed.\n"); 
	my $ttl_base = 0; 
	my $ovl_base = 0; 
	for my $r1 (@$ar) {
		$r1->[1] eq 'CDS' or next; 
		my ($c_chr, $cs, $ce) = (@{$r1->[0]}[0, 3,4]); 
		$ttl_base += ($ce-$cs+1); 
		defined $hr->{$c_chr} or next; 
		for my $r2 (@{$hr->{$c_chr}}) {
			$r2->[1] < $cs and next; 
			$r2->[0] > $ce and last; 
			$ovl_base += ( mathSunhh::min($r2->[1], $ce) - mathSunhh::max($r2->[0], $cs) + 1 ); 
		}
	}
	if ( $ovl_base <= $ttl_base * $ovl_ratio ) {
		return 1; 
	} else {
		return 0; 
	}
}


sub out_gL {
	my $ar = shift;
	for my $r1 (@$ar) {
		print STDOUT join("\t", @{$r1->[0]})."\n";
	}
	return 0;
}
