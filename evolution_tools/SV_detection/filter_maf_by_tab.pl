#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;

!@ARGV and die "perl $0 align2w38_anc_fix.maf.blasttab.good q.fa.kl align2w38_anc_fix.maf > align2w38_anc_fix.good.maf\n";

my $tabF=shift;
my $qKLF=shift;
my %qKL;
for (&fileSunhh::load_tabFile($qKLF)) {
  $_->[0] eq 'key' and next;
  $qKL{$_->[0]} = $_->[1];
}

my %wanted;
for my $l0 (&fileSunhh::load_tabFile($tabF)) {
  my @r = ($l0->[1], $l0->[8]-1, $l0->[9]-$l0->[8]+1, '+');
  if ($l0->[6]-1 == $l0->[7]) {
    my @q = ($l0->[0], $l0->[6]-1, $l0->[7]-$l0->[6]+1, '+');
    my $k=join(" ", @r, @q);
    $wanted{$k} = 1;
    defined $qKL{$l0->[0]} or die "bad 1: |$l0->[0]|\n";
    @q = ($l0->[0], $qKL{$l0->[0]}-$l0->[6], $l0->[6]-$l0->[7]+1 ,'-');
    $k=join(" ", @r, @q);
    $wanted{$k} = 1;
  } elsif ($l0->[6] > $l0->[7]) {
    defined $qKL{$l0->[0]} or die "bad 1: |$l0->[0]|\n";
    my @q = ($l0->[0], $qKL{$l0->[0]}-$l0->[6], $l0->[6]-$l0->[7]+1 ,'-');
    my $k=join(" ", @r, @q);
    $wanted{$k} = 1;
  } elsif ($l0->[6] == $l0->[7]) {
    my @q = ($l0->[0], $l0->[6]-1, $l0->[7]-$l0->[6]+1, '+');
    my $k=join(" ", @r, @q);
    $wanted{$k} = 1;
    defined $qKL{$l0->[0]} or die "bad 1: |$l0->[0]|\n";
    @q = ($l0->[0], $qKL{$l0->[0]}-$l0->[6], $l0->[6]-$l0->[7]+1 ,'-');
    $k=join(" ", @r, @q);
    $wanted{$k} = 1;
  } elsif ($l0->[6]+1 == $l0->[7]) {
    my @q = ($l0->[0], $l0->[6]-1, $l0->[7]-$l0->[6]+1, '+');
    my $k=join(" ", @r, @q);
    $wanted{$k} = 1;
    @q = ($l0->[0], $qKL{$l0->[0]}-$l0->[6], 0, '-');
    $k=join(" ", @r, @q);
    $wanted{$k} = 1;
  } else {
    my @q = ($l0->[0], $l0->[6]-1, $l0->[7]-$l0->[6]+1, '+');
    my $k=join(" ", @r, @q);
    $wanted{$k} = 1;
  }
}

while (<>) {
  chomp;
  if (m!^\s*#!) {
    print STDOUT "$_\n";
    next;
  }
  my $l2=<>; chomp($l2);
  my $l3=<>; chomp($l3);
  my $l4=<>; chomp($l4);
  m!^a\s! or die "bad2:$_\n";
  $l2 =~ m!^s\s+(\S+)\s+(\d+)\s+(\d+)\s+([+-])\s! or die "bad3:$l2\n";
  my @r1=($1,$2,$3,$4);
  $l3 =~ m!^s\s+(\S+)\s+(\d+)\s+(\d+)\s+([+-])\s! or die "bad4:$l3\n";
  my @q1=($1,$2,$3,$4);
  my $k=join(" ", @r1, @q1);
  defined $wanted{$k} or next;
  print "$_\n$l2\n$l3\n$l4\n";
}

