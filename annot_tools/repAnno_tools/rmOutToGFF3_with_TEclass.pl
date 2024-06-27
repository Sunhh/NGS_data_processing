#!/usr/bin/perl
use strict;
use warnings;

if (@ARGV < 1) {
  die "Usage: $0 <RepeatMasker.out file>\n";
}

my $file = $ARGV[0];
open my $in, '<', $file or die "Cannot open $file: $!\n";

while (<$in>) {
  next if /^#/;
  next unless /^\s*\d+/;

  my @fields = split ' ';
  my ($score, $div, $del, $ins, $query, $q_start, $q_end, $q_left, $strand,
    $repeat, $class, $r_start, $r_end, $r_left, $id) = @fields;

  if ($strand eq "C") {
    $strand = "-";
    ($r_start, $r_end) = ($r_left, $r_end);
  } else {
    $strand = "+";
  }

  # my $attributes = "ID=$id;Target=$repeat $r_start $r_end;Class=$class";
  my $attributes = "Target=$repeat $r_start $r_end;Class=$class";
  my $source = "RepeatMasker";
  my $type = "dispersed_repeat";

  print join("\t", $query, $source, $type, $q_start, $q_end, $score, $strand, ".", $attributes), "\n";
}

close $in;
