#!/usr/bin/perl
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;
use mathSunhh;
my $ms_obj = mathSunhh->new();
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "in_blk:s", # ov1/all_type.cds.blk
  "in_grp:s", # ov1/comb.grp2.novl_loc.syn.grp_tbl.filtered
  "max_ovlLen:i", # 0
);

my $htxt = <<HH;
####################################################################################################
# perl $0  -in_grp  ov1/comb.grp2.novl_loc.syn.grp_tbl.filtered  -in_blk ov1/all_type.cds.blk > ov1/comb.grp3
# 
#  -max_ovlLen      [0]
####################################################################################################
HH

for my $k1 (qw/in_grp in_blk/) {
  defined $opts{$k1} or &LogInforSunhh::usage("\n-$k1 is required!\n\n$htxt");
}

$opts{'max_ovlLen'} //= 0;

my ($gene2loc, $chr2gene) = &load_blk($opts{'in_blk'});


{
  my $fh1 = &openFH($opts{'in_grp'}, '<');
  while (<$fh1>) {
    chomp;
    my @ta=split(/\t/, $_);
    my @kept_genes;
    for (my $i=2; $i<@ta; $i++) {
      $ta[$i] =~ m!^\S+:\d+\-\d+:[+-]$! or do { push(@kept_genes, $ta[$i]); next; }; # Only pick a Qloc gene for comparison.
      my $is_bad = 0;
      # Qloc gene.
      for (my $j=2; $j<@ta; $j++) {
        $j == $i and next;
        $ta[$j] =~ m!^\S+:\d+\-\d+:[+-]$! and next; # Only compare this Qloc with predicted genes.
        # predicted gene.
        defined $gene2loc->{$ta[$i]} or die "|$ta[$i]|\n";
        defined $gene2loc->{$ta[$j]} or die "|$ta[$j]|\n";
        $gene2loc->{$ta[$i]}[0] eq $gene2loc->{$ta[$j]}[0] or next;
        $gene2loc->{$ta[$i]}[1] eq $gene2loc->{$ta[$j]}[1] or next;
        $gene2loc->{$ta[$i]}[3][0] > $gene2loc->{$ta[$j]}[3][1] and next;
        $gene2loc->{$ta[$i]}[3][1] < $gene2loc->{$ta[$j]}[3][0] and next;
        my ($ovlLen_nonDup, $ovlCnt_mayDup, $ovlLocAR) = $ms_obj->compare_number_list( $gene2loc->{$ta[$i]}[2], $gene2loc->{$ta[$j]}[2], 'compare' => 'ovl', 'sort'=>'0' );
        if ($ovlLen_nonDup > $opts{'max_ovlLen'}) {
          # I should remove $i.
          $is_bad = 1;
          last;
        }
      }
      $is_bad == 1 and next;
      push(@kept_genes, $ta[$i]);
    }
    $ta[1] = scalar(@kept_genes);
    print STDOUT join("\t", $ta[0], $ta[1], @kept_genes)."\n";
  }
  close($fh1);
}
# GrpSyn_000001   31      CApan:CaUC09G173540.1   CApan:CaUC09G173550.1   CApan:CaUC09G173640.1                   CApan:CaUC09G173650.1                   CApan:CaUC09G173660.1                   CApan:CaUC09G173670.1
# GrpSyn_000002   21      CApan:CaUC03G061460.1   CApan:CaUC03G061520.1   CApan:CaUC03G061550.1                   CApan:CaUC03G061560.1                   CApan:Ciama_Chr03:24544226-24545011:+   CApan:Ciama_Chr03:30256016-30256508:+
# GrpSyn_000003   19      CApan:CaUC10G186420.1   CApan:CaUC10G186430.1   CApan:Ciama_Chr10:12115704-12116238:+   CLpan:Cla97C10G191540.1                 CLpan:Cla97C10G191550.1                 CLpan:Cla97C10G191553.1
# GrpSyn_000004   19      CCpan:CcUC04G060480.1   CCpan:CcUC04G060490.1   CCpan:CicolChr04:3296818-3299346:-      CCpan:CicolChr04:7989701-7990121:-      CCpan:CicolChr04:8006641-8007081:-      CLpan:Cla97C04G069310.1

# Return: (\%gene2loc, \%chr2gene);
#   $gene2loc{mrnaID}     = [ chrID, chrStrand, [ [s1,e1], [s2,e2], ... ], [s, e] ] # small to large; 
#   $chr2gene{chrID}{str} = [ [start, end, mrnaID], [], ... ]
sub load_blk {
  my ($fn) = @_; # $fn is already sorted, with the first exon at larger position for minus strand genes.
  my %back;
  my %b_c2g;
  for my $l1 (&fileSunhh::load_tabFile($fn)) {
    $back{$l1->[0]} = [ $l1->[1], $l1->[2], [], []]; # $back{$mrnaID} = [ chrID, chrStrand, [ [s1,e1], [s2,e2], ... ], [s, e] ]
    my ($s, $e);
    for my $a1 (split(/;/, $l1->[3])) {
      $a1 =~ m!^(\d+),(\d+)$! or &stopErr("[Err] Bad loc [$a1]\n");
      push(@{$back{$l1->[0]}[2]}, [$1, $2]);
      $s //= $1; $s > $1 and $s = $1;
      $e //= $2; $e < $2 and $e = $2;
    }
    $back{$l1->[0]}[3] = [$s, $e];
    @{$back{$l1->[0]}[2]} = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$back{$l1->[0]}[2]};
    push(@{$b_c2g{$l1->[1]}{$l1->[2]}}, [$s, $e, $l1->[0]]); # $b_c2g{$chrID}{$strand} = [ [ start, end, mrnaID ], [], ... ];
  }
  for my $chrID (keys %b_c2g) {
    for my $str (keys %{$b_c2g{$chrID}}) {
      @{$b_c2g{$chrID}{$str}} = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$b_c2g{$chrID}{$str}};
    }
  }
  return(\%back, \%b_c2g);
}# load_blk()
# Sunhh@noname:/data/Sunhh/wmpan1_ortho_genes/by_liftoff$ head -4 ov1/all_type.cds.blk | less -S
# CApan:CaUC00G218630.1   Ciama_Chr00     +       233297,234241;234551,234643
# CApan:CaUC00G218640.1   Ciama_Chr00     -       251253,253637
# CApan:CaUC00G218650.1   Ciama_Chr00     +       362928,362933;363002,363021;363232,364222;365053,365670



