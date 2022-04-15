#!/usr/bin/perl
# [3/28/2022] CDS blocks are required for further use.
use strict;
use warnings;

!@ARGV and die "perl $0 sp_tag in.psl > in.psl.gene_loc\n";

my $tag = shift;

my $min_cov   = 0.9;
my $min_ident = 0.95;

my @a1;
while (<>) {
  chomp;
  m!^\d+! or next;
  my @ta=split(/\t/, $_);
  $ta[9] =~ m!^$tag:! and next;
  $ta[12]-$ta[11] >= $min_cov * $ta[10] or next;
  $ta[0] >= $min_ident * ($ta[0] + $ta[1]) or next;
  my ($qblks, $tblks) = &psl_line_to_Blks(\@ta);
  push(@a1, [$ta[9], "$tag:$ta[13]:".($ta[15]+1)."-$ta[16]:$ta[8]", $ta[0], $ta[12]-$ta[11], $ta[13], $ta[8], $tblks]);
}
@a1 = sort { $b->[2] <=> $a->[2] || $a->[3] <=> $b->[3] } @a1;
my %h1;
for (@a1) {
  defined $h1{$_->[0]} and next;
  $h1{$_->[0]} = 1;
  print STDOUT join("\t", @$_)."\n";
}

sub psl_line_to_Blks {
  my ($ar) = @_;
  my @blockSizes = split(/,/, $ar->[18]);
  my @qStarts    = split(/,/, $ar->[19]);
  my @tStarts    = split(/,/, $ar->[20]);
  my (@qblks, @tblks);
  if      ($ar->[8] eq '+') {
    for (my $i=0; $i<@qStarts; $i++) {
      push(@qblks, join(",", $qStarts[$i]+1, $qStarts[$i]+$blockSizes[$i]));
      push(@tblks, join(",", $tStarts[$i]+1, $tStarts[$i]+$blockSizes[$i]));
    }
  } elsif ($ar->[8] eq '-') {
    for (my $i=0; $i<@qStarts; $i++) {
      push( @qblks, join(",", $ar->[10]-$qStarts[$i]-$blockSizes[$i]+1, $ar->[10]-$qStarts[$i]) );
      push( @tblks, join(",", $tStarts[$i]+1, $tStarts[$i]+$blockSizes[$i]) );
    }
    @qblks = reverse(@qblks);
    @tblks = reverse(@tblks);
  } else {
    die "bad str [$ar->[8]] in line: @$ar\n";
  }
  return(join(";", @qblks), join(";", @tblks));
}# psl_line_to_TBlks()

