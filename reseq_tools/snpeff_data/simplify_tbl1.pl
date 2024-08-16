#!/usr/bin/perl
# 7/11/2024: v2. Add comments. Instead of a merge line for each SV, output a melt table with only one SV vs. one gene in one line.
use strict;
use warnings;
use fileSunhh;

!@ARGV and die "perl $0 map.geneID_mrnaID in_snpeff.tbl1 > out\n";

my $fn_mapID = shift; # /data/Sunhh/src/evolution/snpeff/snpEff/data/97pan/map.geneID_mrnaID
my %g2m;
{
#warn "[Msg] Loading mapping file [$fn_mapID]\n";
for my $l (&fileSunhh::load_tabFile($fn_mapID)) { $g2m{$l->[0]} = $l->[1]; }
#warn "[Msg] Loaded\n";
}

my $fn_svClass = '/home/Sunhh/tools/github/NGS_data_processing/reseq_tools/snpeff_data/simple_sv_class';
my $gen_cn = 6; # 'ANN[*].FEATUREID';
my $eff_cn = 7; # 'ANN[*].EFFECT';
my $dis_cn = 8; # '[*].DISTANCE';

# 'allInDel-ann_ud2k.tbl1' / 'in_snpeff.tbl1'; Output from command:
## snpsift extractFields -s ',' -e '.' - \
##   CHROM POS ID SVLEN ALT "ANN[*].ALLELE" "ANN[*].FEATUREID" "ANN[*].EFFECT" "ANN[*].DISTANCE" \
##   > allInDel-ann_ud2k.tbl1 ;

my %svtype; # {snpeff_type}{simple_1/order_1/...} = value_text;
{
  my @d1 = &fileSunhh::load_tabFile($fn_svClass);
  my @h1 = @{shift(@d1)};;
  for my $l (@d1) {
    for (my $i=1;$i<@h1;$i++) {
      $svtype{$l->[0]}{$h1[$i]} = $l->[$i];
    }
  }
}


while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  if ($ta[0] eq 'CHROM') {
    print STDOUT join("\t", @ta[0,1,2,3], qw/gene jnEffect distance simple_1 oriEffect/)."\n";
    next;
  }
  my @g1 = split(/,/, $ta[$gen_cn]); # gene names.
  for my $g1t (@g1) { defined $g2m{$g1t} and $g1t=$g2m{$g1t}; }
  my @e1 = split(/,/, $ta[$eff_cn]); # effects.
  my @d1 = split(/,/, $ta[$dis_cn]); # distances.
  my %eff1;
  ### {geneID} = [
  ###   [ simple_1_eff, simple_2_eff, [raw_eff_1, raw_eff_2, ...
  ###                                 ], order_1_rank, distance
  ###   ],
  ###   [], ...
  ### ];
  # For each sv-gene_effect pair, keep only the effects with the highest 'order_1' rank;
  #   Group them by genes.
  for (my $i=0; $i<@g1; $i++) {
    $g1[$i] eq '.' and next; # The feature ID is '.' for 'feature_ablation' because it uses gene name.
    my @eff2;
    # Split effects like 'frameshift_variant&splice_region_variant' to keep only one effect.
    for my $e2 (split(/\&/, $e1[$i])) {
      $e2 eq '.' and next;
      defined $svtype{$e2} or die "[Err] Failed to rank effect [E2|$ta[2]|$e2]\n";
      $svtype{$e2}{'simple_1'} eq 'intergenic' and next; # Skip intergenic region type.
      if (scalar(@eff2) == 0 or $svtype{$e2}{'order_1'} < $eff2[3]) {
        @eff2 = ($svtype{$e2}{'simple_1'}, $svtype{$e2}{'simple_2'}, [$e2], $svtype{$e2}{'order_1'}, $svtype{$e2}{'order_3'});
      } elsif ($svtype{$e2}{'order_1'} == $eff2[3]) {
        push(@{$eff2[2]}, $e2);
      }
    }
    scalar(@eff2) > 0 and push(@{$eff1{$g1[$i]}}, [@eff2, $d1[$i]]);
  }
  # For each sv-gene pair, keep only one effect with highest 'order_3' rank;
  #   Select the only one of the detailed effect for each gene.
  #   Remove intergenic type.
  #   Remove duplicated effects within the same gene.
  my @eff4;
  for my $tg1 (keys %eff1) {
    @{$eff1{$tg1}} = sort {$a->[3] <=> $b->[3]} @{$eff1{$tg1}}; # Sort effects from the same gene.
    my @tef1 = @{$eff1{$tg1}[0]};
    my $tord3 = $tef1[4]; # 'order_3';
    if (scalar(@eff4) == 0) {
      @eff4 = ([join(";", $tg1, $tef1[0], $tef1[2][0], $tef1[5]), $tord3, $tg1, $tef1[5], $tef1[0], $tef1[2][0]]);
    } elsif ($tord3 < $eff4[-1][1]) {
      @eff4 = ([join(";", $tg1, $tef1[0], $tef1[2][0], $tef1[5]), $tord3, $tg1, $tef1[5], $tef1[0], $tef1[2][0]]);
    } elsif ($tord3 == $eff4[-1][1]) {
      push(@eff4, [join(";", $tg1, $tef1[0], $tef1[2][0], $tef1[5]), $tord3, $tg1, $tef1[5], $tef1[0], $tef1[2][0]]);
    }
  }
  if (scalar(@eff4) > 0) {
    for my $t2 (sort {$a->[0] cmp $b->[0]} @eff4) {
      print join("\t", @ta[0,1,2,3], $t2->[2], $t2->[0], $t2->[3], $t2->[4], $t2->[5])."\n";
    }
  } else {
    print join("\t", @ta[0,1,2,3], qw/. . . . ./)."\n";
  }
}

