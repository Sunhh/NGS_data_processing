#!/usr/bin/perl
use strict; 
use warnings; 

use LogInforSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"debug!", 
); 

-t and !@ARGV and die "perl $0 dup_aug1_aug2.gff3 > nonDup_aug1_aug2.gff3\n"; 

my %have_model; 

# ##gff-version 3
# S401083_pilon   AUGUSTUS        match   113     6130    0.39    -       .       ID=0:g1.t1;
# S401083_pilon   AUGUSTUS        match_part      113     472     1       -       0       ID=0:g1.t1.c1;Parent=0:g1.t1;
# S401083_pilon   AUGUSTUS        match_part      595     807     0.99    -       0       ID=0:g1.t1.c2;Parent=0:g1.t1;
# S401083_pilon   AUGUSTUS        match_part      1269    1412    0.91    -       0       ID=0:g1.t1.c3;Parent=0:g1.t1;

my @geneLines; 
while (<>) {
	m/^#/ and next; 
	m/^\s*$/ and next; 
	m/^>/ and last; 
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($ta[2] eq 'match') {
		if (scalar(@geneLines) > 0) {
			&is_have(\@geneLines, \%have_model) or &out_gL(\@geneLines); 
		}
		@geneLines = (); 
		push(@geneLines, [[@ta]]); 
	} elsif ( $ta[2] eq 'match_part' ) {
		push(@geneLines, [[@ta]]); 
	} else {
		die; 
	}
}

sub is_have {
	my ($ar, $hr) = @_; 
	my $scfID = $ar->[0][0][0]; 
	my $is_have = 1; 
	my ($gS, $gE, $gStr) = @{$ar->[0][0]}[3,4,6]; 
	my $gKey = join("\t", $gS, $gE, $gStr); 
	my $had_ta8 = ''; 
	
	if ( defined $hr->{$scfID} ) {
		if ( defined $hr->{$scfID}{$gKey} ) {
			$is_have = 0; 
			CHK_MODEL: 
			for my $tr (@{$hr->{$scfID}{$gKey}}) {
				my $tmp_have = 1; 
				CHK_EXON: 
				for (my $i=0; $i<@{$tr}; $i++) {
					$ar->[$i][0][3] eq $tr->[$i][0][3] or do { $tmp_have = 0; last CHK_EXON; } ; 
					$ar->[$i][0][4] eq $tr->[$i][0][4] or do { $tmp_have = 0; last CHK_EXON; } ; 
					$ar->[$i][0][6] eq $tr->[$i][0][6] or do { $tmp_have = 0; last CHK_EXON; } ; 
				}
				$tmp_have == 1 and do { $is_have = 1; $had_ta8 = $tr->[0][0][8]; last CHK_MODEL; }; 
			}
		} else {
			$is_have = 0; 
		}
	} else {
		$is_have = 0; 
	}
	if ($is_have == 0) {
		push(@{$hr->{$scfID}{$gKey}}, [@$ar]); 
	} else {
		$opts{'debug'} and &tsmsg("[Msg] Skip duplicted gKey [$gKey] for [$ar->[0][0][8]] by [$had_ta8]\n"); 
	}
	
	return $is_have; 
}

if (scalar(@geneLines) > 0) {
	&is_have(\@geneLines, \%have_model) or &out_gL(\@geneLines); 
	@geneLines = (); 
}

sub out_gL {
	my $ar = shift;
	for my $r1 (@$ar) {
		print STDOUT join("\t", @{$r1->[0]})."\n";
	}
	return 0;
}
