#!/usr/bin/perl
# [12/21/2022] Fix out.maf output from AnchorWave.
use strict;
use warnings;

-t and !@ARGV and die "perl $0 AnchorWave_out_align1.maf > fixed.maf\n";

while (<>) {
  m!^#! and do { print; next; };
  m!^a\s! or die "[Err] Unexpected format of MAF:$_\n";
  my $s1=<>;
  my $s2=<>;
  my $blank=<>;
  $blank =~ m!^$! or die "[Err] Unexpected blank line: $blank\n";
  $s1 = &fix_strPos($s1);
  $s2 = &fix_strPos($s2);
  print "$_$s1$s2$blank";
}

sub fix_strPos {
  my ($ss) = @_;
  $ss =~ m!^s\s+\S+\s+(\d+)\s+(\d+)\s+([+-])\s+(\d+)! or die "[Err] Unexpected s line: $ss\n";
  if ($3 eq '-') {
    my $newS = $4-($1+$2);
    $ss =~ s!^(s\s+\S+\s+)\d+(\s+)!$1$newS$2! or die "[Err] Failed to change: $ss\n";
  } elsif ($3 eq '+') {
    ;
  } else {
    die "[Err] Impossible strand $3: $ss\n";
  }
  return($ss);
}# fix_strPos()

