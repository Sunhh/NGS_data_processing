#!/usr/bin/perl -w
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"flank_extend:i", # [0]
	"inAnnotLis:s", 
	"inLocLis:s", 
	"tagLoc_colN:i", # -9999
); 

$opts{'tagLoc_colN'} //= -9999; 

sub usage {
	print <<HH; 
################################################################################
# perl $0 -inLocLis apple_noN_w10ks10k.p1.compare.grp_2to1.slct -inAnnotLis Malus_x_domestica.v1.0-primary.transcripts.gff3.annot
# 
# -help 
# -inLocLis        Required. File with location selected. 
#                    Format: chr10   30320001        30330000        10000   10000
#                            chrID   chrS            chrE            lenSpan BpCnt
#   -tagLoc_colN   [$opts{'tagLoc_colN'}] -9999 means not assigned. 
# -inAnnotLis      Required. File with annotation of genes. 
#                    Format: chr10   73429   74098   -          MDP0000184598   "GO:0031072"    predicted with fgenesh_arab   Lg10_Scaffold1:73429..74098-
#                            chrID   chrS    chrE    chrStrand       ProtID     GO              Note                          Scaffold
# -flank_extend    [0] Length of flanking extensions for the inLocLis. 
################################################################################
HH
	exit 1; 
}
$opts{'help'} and &usage(); 
( defined $opts{'inLocLis'} and defined $opts{'inAnnotLis'} ) or &usage(); 
$opts{'flank_extend'} //= 0; 

my $inLocFh = &openFH( $opts{'inLocLis'}, '<' ); 
my $inAnnFh = &openFH( $opts{'inAnnotLis'}, '<' ); 
my %geneInChr; 
my %chrLoc; 

########### Read in the loc file. 

my %loc2tag; 

while (<$inLocFh>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($. == 1 and $ta[0] =~ m/^chrID$/i) {
		next; 
	}
	my ($chrID, $chrS, $chrE) = @ta; 
	push(@{$chrLoc{$chrID}}, [$chrS, $chrE]); 
	my $loc_ID = "$chrID:$chrS-$chrE"; 
	if ($opts{'tagLoc_colN'} != -9999) {
		push(@{$loc2tag{$loc_ID}}, $ta[$opts{'tagLoc_colN'}]); 
	} else {
		push(@{$loc2tag{$loc_ID}}, $loc_ID); 
	}
} 
for my $chrID (sort keys %chrLoc) {
	@{$chrLoc{$chrID}} = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$chrLoc{$chrID}} ; 
}

while (<$inAnnFh>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($. == 1 and $ta[0] =~ m/^chrID$/i) {
		print STDOUT join("\t", $_, 'Ovl_Loc')."\n"; 
		next; 
	}
	my ($chrID, $chrS, $chrE, $chrStr) = @ta; 
	my $is_in  = 0; 
	my @ovl_locs; 
	for my $se_ar ( @{$chrLoc{$chrID}} ) {
		my ($s, $e) = @$se_ar; 
		$s - $opts{'flank_extend'} > $chrE and last; 
		$e + $opts{'flank_extend'} < $chrS and next; 
		$is_in = 1; 
		my $loc_ID = "$chrID:$s-$e"; 
		defined $loc2tag{$loc_ID} or die "|$loc_ID|\n"; 
		push(@ovl_locs, @{$loc2tag{$loc_ID}}); 
	}
	if ($is_in == 1) {
		print STDOUT join("\t", $_, join(";", @ovl_locs))."\n"; 
	}
}
