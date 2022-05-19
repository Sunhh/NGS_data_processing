#!/usr/bin/perl
use strict;
use warnings;

-t and !@ARGV and die "perl $0 mcscan_run_p/ASM1002_97103_vs_ASM1003_PI537277.tbl.orth.gene_pairs > mcscan_run_p/ASM1002_97103_vs_ASM1003_PI537277.tbl.orth.gene_pairs.OG\n";

my %h;
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'BlkID' and next;
  &update_best(\%h, $ta[1], $ta[2], $ta[3]);
  &update_best(\%h, $ta[2], $ta[1], $ta[3]);
}

my %has;
for my $id1 (sort {$h{$a}[1] <=> $h{$b}[1] || scalar(@{$h{$b}[0]}) <=> scalar(@{$h{$b}[0]}) || $a cmp $b} keys %h) {
  defined $has{$id1} and next;
  $has{$id1} = 1;
  # @{$h{$id1}[0]}
  my @good_id2;
  for my $id2 (@{$h{$id1}[0]}) {
    defined $has{$id2} and next;
    push(@good_id2, $id2);
    $has{$id2} = 1;
  }
  scalar(@good_id2) > 0 or next;
  print join("\t", $id1, join(";;", @good_id2), $h{$id1}[1])."\n";
}

sub update_best {
  my ($hR, $v1, $v2, $v3) = @_;
  if (defined $hR->{$v1}) {
    if      ( $hR->{$v1}[1] > $v3 ) {
      $hR->{$v1} = [ [$v2], $v3 ];
    } elsif ( $hR->{$v1}[1] == $v3 ) {
      push(@{$hR->{$v1}[0]}, $v2);
    } else {
      # leave hR unchanged.
    }
  } else {
    $hR->{$v1} = [ [$v2], $v3 ];
  }
  return;
}# update_best()

# BlkID	Gene1	Gene2	Ks
# 29	CcUC01G000020.1	Cla97C01G000010.2	0.0464
# 29	CcUC01G000030.1	Cla97C01G000020.1	0.0275
# 29	CcUC01G000040.1	Cla97C01G000030.2	0.0651
