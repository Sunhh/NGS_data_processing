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
  $ta[1] eq 'Final' and do { print "$_\n"; next; };
  $ta[11] < 1 and next;
  if ($ta[11] == 1) {
    $ta[1] eq $expect_Final and $ta[1] eq $expect_Expan and print "$_\n";
  } else {
    $ta[1] eq $expect_Final and print "$_\n";
  }
}

