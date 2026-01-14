#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0 SP in.gff3.JnLoc > out_mcscanx.gff\n";

my $spName = shift;
length($spName) == 2 or die "Spedies prefix should be two characters.";

my @o1;
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'mrnaID' and next;
  push(@o1, [@ta[2,0], $ta[6]-1, $ta[7], $ta[2], $ta[5]]);
}
@o1 = sort {$a->[0] cmp $b->[0] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3]} @o1;
my %chr2num;
my $chrcnt = 0;
for (@o1) {
  # $_->[4] = $_->[0];
  unless (defined $chr2num{$_->[0]}) {
    $chrcnt ++;
    $chr2num{$_->[0]} = $spName . $chrcnt;
  }
  $_->[0] = $chr2num{$_->[0]};
  print STDOUT join("\t", @$_)."\n";
}

