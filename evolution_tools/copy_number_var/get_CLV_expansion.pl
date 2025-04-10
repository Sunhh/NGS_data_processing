#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my %opts;
GetOptions(\%opts,
  "final_label:s",
  "expan_label:s",
  "help!",
);

$opts{'final_label'} //= 'CLV_high';
$opts{'expan_label'} //= 'CLV_high';

my $htxt = <<HH;
perl $0 CA_to_CLV.tbl [-final_label CLV_high  -expan_label CLV_high]
HH

-t and !@ARGV and die $htxt;

my $expect_Final = $opts{'final_label'};
my $expect_Expan = $opts{'expan_label'};

while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[1] eq 'Final' and do { print join("\t", qw/Gene base_size Mode_all Final Expansion Contraction/)."\n"; next; };
  $ta[11] < 1 and next;
  if ($ta[11] == 1) {
    if ($ta[1] eq $expect_Final and $ta[2] eq $expect_Expan) {
      # >= 2 gene family size is enriched in 'expect_Final';
      print join("\t", $ta[0], 1, @ta[11,1,2,3])."\n";
    }
  } else {
    if ($ta[1] eq $expect_Final and $ta[2] eq $expect_Expan) {
      # > Mode_all ($ta[11]) is enriched in 'expect_Final'
      print join("\t", $ta[0], $ta[11], @ta[11,1,2,3])."\n";
    } elsif ($ta[1] eq $expect_Final) {
      # < Mode_all ($ta[11]) contraction is enriched in the other group.
      print join("\t", $ta[0], $ta[11]-1, @ta[11,1,2,3])."\n";
    }
  }
}

