#!/usr/bin/perl
# [3/17/2022] Any CDS region exceeding the truncated novel contig part will be trimmed. The novel contig region should be continuous on the ext-contig.
# [5/12/2022] If a gene has two CDS on two different novel contigs belonging to the same extended (original) contig, I would choose the one with longest CDS region.
use strict;
use warnings;

!@ARGV and die "perl $0 CL.ext2nov.tbl in_ext.gff3 > trimmed_ext.gff3\n";

my $f1_tbl = shift;
my $f2_gff = shift;

open F1,'<',"$f1_tbl" or die;
# [Sunhh@panda rmTE]$ head -4 /data/Sunhh/watermelon/08.pan_phase_1/CL/04_novel_seqs/01_gaolei_approach/final/CL.ext2nov.tbl
# 182M__NODE_54395__1_576_ext     1       576     +       182M_NODE_54395F        182M
# 182M__NODE_55106__1_505_ext     1       505     +       182M_NODE_55106F        182M
# 213__NODE_18317__1_4868_ext     1       525     +       213_NODE_18317__1-525H  213
my %ext2nov;
while (<F1>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[3] eq "+" or die "I accept only + strand only.\n";
  push(@{$ext2nov{$ta[0]}}, [@ta[1,2], $ta[4]]);
}
close F1;
for (keys %ext2nov) {
  @{$ext2nov{$_}} = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$ext2nov{$_}};
}

open F2,'<',"$f2_gff" or die;
# [Sunhh@panda rmTE]$ head -4 final_cln.gff3
# Zz61__NODE_26173__1_2578_ext    maker   gene    2303    2533    .       -       .       ID=maker-s032383-pred_gff_augustus_masked-gene-0.0;Name=maker-s032383-pred_gff_augustus_masked-gene-0.0;score=0.51;target_length=2578
# Zz61__NODE_26173__1_2578_ext    maker   mRNA    2303    2533    .       -       .       ID=maker-s032383-pred_gff_augustus_masked-gene-0.0-mRNA-1;Parent=maker-s032383-pred_gff_augustus_masked-gene-0.0;Name=maker-s032383
# Zz61__NODE_26173__1_2578_ext    maker   exon    2303    2533    .       -       .       ID=maker-s032383-pred_gff_augustus_masked-gene-0.0-mRNA-1:1;Parent=maker-s032383-pred_gff_augustus_masked-gene-0.0-mRNA-1
# Zz61__NODE_26173__1_2578_ext    maker   CDS     2303    2533    .       -       0       ID=maker-s032383-pred_gff_augustus_masked-gene-0.0-mRNA-1:cds;Parent=maker-s032383-pred_gff_augustus_masked-gene-0.0-mRNA-1
my %info;
while (<F2>) {
  m!^\s*(#|$)! and do {print; next;};
  chomp;
  my @ta=split(/\t/, $_);
  if (!defined $ext2nov{$ta[0]}) {
    warn "[Wrn] Skipping line: $_\n";
    print STDOUT join("\t", @ta)."\n";
    next;
  }
  my ($oldS, $oldE) = @ta[3,4];
  my @newSE;
  for my $a1 (@{$ext2nov{$ta[0]}}) {
    $a1->[0] > $oldE and last;
    $a1->[1] < $oldS and next;
    my ($newS, $newE);
    if ($a1->[0] <= $oldS) {
      $newS //= $oldS;
    } else {
      $newS //= $a1->[0];
    }
    if ($a1->[1] >= $oldE) {
      $newE //= $oldE;
    } else {
      $newE //= $a1->[1];
    }
    $newE >= $newS or die "why here 1!\n";
    push(@newSE, [$newS, $newE, $a1->[2]]);
    next;
  }
  scalar(@newSE) > 0 or do { warn "[Wrn] Skip element out of range: $_\n"; next; };
  for (my $i=0; $i<@newSE; $i++) {
    my ($newS, $newE, $newCtgID) = @{$newSE[$i]};
    @ta[3,4] = ($newS, $newE);
    if      ($ta[2] eq "gene") {
      $ta[8] =~ m!ID=([^\s;]+)! or die "gene:$ta[8]\n";
      $info{'gene'}{$1} = [[@ta]];
    } elsif ($ta[2] eq "mRNA") {
      $ta[8] =~ m!ID=([^\s;]+)! or die "mrnaID:$ta[8]\n";
      my $mid = $1;
      $ta[8] =~ m!Parent=([^\s;]+)! or die "mrnaParent:$ta[8]\n";
      my $gid = $1;
      $info{'mrna'}{$mid} = [[@ta], $gid, $.];
    } elsif ($ta[2] eq "exon") {
      $ta[8] =~ m!Parent=([^\s;]+)! or die "exonParent:$ta[8]\n";
      my $mid = $1;
      push(@{$info{'exon'}{$mid}}, [@ta]);
      &update_SE(\%info, $mid, $newS, $newE);
    } elsif ($ta[2] eq "CDS") {
      $ta[8] =~ m!Parent=([^\s;]+)! or die "cdsParent:$ta[8]\n";
      my $mid = $1;
      push(@{$info{'cds'}{$mid}}, [@ta]);
      &update_SE(\%info, $mid, $newS, $newE);
    } elsif ($ta[2] eq "five_prime_UTR") {
      $ta[8] =~ m!Parent=([^\s;]+)! or die "5utrParent:$ta[8]\n";
      my $mid = $1;
      push(@{$info{'5utr'}{$mid}}, [@ta]);
      &update_SE(\%info, $mid, $newS, $newE);
    } elsif ($ta[2] eq "three_prime_UTR") {
      $ta[8] =~ m!Parent=([^\s;]+)! or die "3utrParent:$ta[8]\n";
      my $mid = $1;
      push(@{$info{'3utr'}{$mid}}, [@ta]);
      &update_SE(\%info, $mid, $newS, $newE);
    } else {
      die "[Err] Bad line type: $_\n";
    }
  }
}
close F2;
my %eleRank = qw(
gene    1
mrna    2
exon    3
five_prime_utr  4
cds             5
three_prime_utr 6
);

for my $mid (sort { $info{'mrna'}{$a}[2] <=> $info{'mrna'}{$b}[2] } keys %{$info{'mrna'}}) {
  my $gid = $info{'mrna'}{$mid}[1];
  my $str = $info{'mrna'}{$mid}[0][6]; # + / - ;
  my (@o0, @o1);
  defined $info{'SE'}{$mid} or die "[Err] No SE found for [$mid]\n";
  my ($s, $e) = @{$info{'SE'}{$mid}};

  # Check if $s and $e are in the same novel contig and then update ($s, $e) to make them in the same novel contig.
  my $chrID = $info{'mrna'}{$mid}[0][0];
  defined $ext2nov{$chrID} or die "why here 2!\n";
  my $trimID = '';
  my %trimCov; # Count the covered length of child elements of mRNA for novel contig selection.
  for my $a1 (@{$ext2nov{$chrID}}) {
    $a1->[0] > $e and last;
    $a1->[1] < $s and next;
    my $blkS = ($s > $a1->[0]) ? $s : $a1->[0] ;
    my $blkE = ($e < $a1->[1]) ? $e : $a1->[1] ;
    my $tk = "$blkS\t$blkE";
    if ($s >= $a1->[0] and $e <= $a1->[1]) {
      # This novel contig contains the whoe gene.
      $trimID = $tk;
      last;
    } else {
      # This gene spans two different novel contigs.
      if (defined $info{'cds'}{$mid}) {
        for my $l1 (@{$info{'cds'}{$mid}}) {
          my $ts = ($a1->[0] > $l1->[3]) ? $a1->[0] : $l1->[3] ;
          my $te = ($a1->[1] < $l1->[4]) ? $a1->[1] : $l1->[4] ;
          $trimCov{$tk} += ($te-$ts+1);
        }
      } else {
        for my $ttype (qw/exon five_prime_utr three_prime_utr/) {
          defined $info{$ttype}{$mid} or next;
          for my $l1 (@{$info{$ttype}{$mid}}) {
            my $ts = ($a1->[0] > $l1->[3]) ? $a1->[0] : $l1->[3] ;
            my $te = ($a1->[1] < $l1->[4]) ? $a1->[1] : $l1->[4] ;
            $trimCov{$tk} += ($te-$ts+1);
          }
        }
      }
    }
  }
  $trimID eq '' and ($trimID) = sort { $trimCov{$b} <=> $trimCov{$a} } keys %trimCov; # The $s and $e are in different novel contigs.
  ($s, $e) = split(/\t/, $trimID);

  # Update the position of each element for output.
  if (defined $info{'gene'}{$gid}) {
    $info{'gene'}{$gid}[0][3] = $s;
    $info{'gene'}{$gid}[0][4] = $e;
    push(@o0, [@{$info{'gene'}{$gid}[0]}]);
  }
  $info{'mrna'}{$mid}[0][3] = $s;
  $info{'mrna'}{$mid}[0][4] = $e;
  push(@o0, [@{$info{'mrna'}{$mid}[0]}]);
  for my $ttype (qw/exon five_prime_utr cds three_prime_utr/) {
    defined $info{$ttype}{$mid} or next;
    for my $l1 (@{$info{$ttype}{$mid}}) {
      $l1->[3] < $s and $l1->[4] >= $s and die "why here 3!\n";
      $l1->[3] <= $e and $l1->[4] > $e and die "why here 4!\n";
      ($l1->[3] >= $s and $l1->[4] <= $e) or next;
      push(@o1, [@$l1]);
    }
  }
  if ($str eq "+") {
    @o1 = sort { $a->[3] <=> $b->[3] || $eleRank{lc($a->[2])} <=> $eleRank{lc($b->[2])} } @o1;
  } elsif ($str eq "-") {
    @o1 = sort { $b->[4] <=> $a->[4] || $eleRank{lc($a->[2])} <=> $eleRank{lc($b->[2])} } @o1;
  } else {
    die "[Err] Unknown str [$str]\n";
  }
  for my $l1 (@o0, @o1) {
    print STDOUT join("\t", @$l1)."\n";
  }
}

sub update_SE {
  my ($hR, $mid, $s, $e) = @_;
  $hR->{'SE'}{$mid} //= [$s, $e];
  $hR->{'SE'}{$mid}[0] > $s and $hR->{'SE'}{$mid}[0] = $s;
  $hR->{'SE'}{$mid}[1] < $e and $hR->{'SE'}{$mid}[1] = $e;
  return;
}# update_SE()

# Zz61__NODE_26173__1_2578_ext    maker   gene    2303    2533    .       -       .       ID=maker-s032383-pred_gff_augustus_masked-gene-0.0;Name=maker-s032383-pred_gff_augustus_masked-gene-0.0;score=0.51;target_length=2578
# Zz61__NODE_26173__1_2578_ext    maker   mRNA    2303    2533    .       -       .       ID=maker-s032383-pred_gff_augustus_masked-gene-0.0-mRNA-1;Parent=maker-s032383-pred_gff_augustus_masked-gene-0.0;Name=maker-s032383
# Zz61__NODE_26173__1_2578_ext    maker   exon    2303    2533    .       -       .       ID=maker-s032383-pred_gff_augustus_masked-gene-0.0-mRNA-1:1;Parent=maker-s032383-pred_gff_augustus_masked-gene-0.0-mRNA-1

