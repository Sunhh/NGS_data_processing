#!/usr/bin/perl
# [3/28/2022] I require the gene IDs in different pan-genome are different.
use strict;
use warnings;

-t and !@ARGV and die "perl $0 output/mapCDS.C?pan.to.C?pan.tbl.gene_map > comb.grp1\n";


my %g2g;

while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'ori_R_ID' and next;
  my ($fromID, $toID) = @ta[0,1];
  $g2g{$fromID}{$toID} = 1;
  $g2g{$toID}{$fromID} = 1;
}

# Group groups.
my @grp1 = @{ &grp_g2g(\%g2g) };
for (@grp1) { @$_ = sort @$_; }

# Output
my $nLen = length(scalar(@grp1));
$nLen ++;
my $gCnt = 0;
for my $tg (sort { scalar(@$b) <=> scalar(@$a) || $a->[0] cmp $b->[0] } @grp1) {
  $gCnt ++;
  # my $grpSize = scalar(@$tg);
  # print STDOUT join("\t", $grpSize, @$tg)."\n";
  my $grpID = sprintf("Grp_%0${nLen}d", $gCnt);
  print STDOUT join("\t", $grpID, scalar(@$tg), @$tg)."\n";
}

sub grp_g2g {
  my ($hr) = @_;
  my %done_g;
  my @back_grp;
  for my $g1 (keys %$hr) {
    defined $done_g{$g1} and next;
    my @g2 = &sub_g2g($hr, [$g1]);
    push(@back_grp, [ $g1, @g2 ]);
    $done_g{$g1} = 1;
    for my $tg2 (@g2) {
      $done_g{$tg2} = 1;
    }
  }
  return(\@back_grp);
}# grp_g2g()
sub sub_g2g {
  my ($hr, $has_ar) = @_;
  my %in_g;
  for my $tg (@$has_ar) { $in_g{$tg} = 1; }
  my %new_g;
  for my $has_g (@$has_ar) {
    for my $tg (keys %{$hr->{$has_g}}) {
      defined $in_g{$tg} and next;
      defined $new_g{$tg} and next;
      $new_g{$tg} = 1;
    }
  }
  my @new_g_arr = keys %new_g;
  scalar(@new_g_arr) == 0 and return();
  my @g3 = &sub_g2g($hr, [@$has_ar, @new_g_arr]);
  push(@new_g_arr, @g3);
  return(@new_g_arr);
}# sub_g2g()

