#!/usr/bin/perl
use strict; 
use warnings;

-t and !@ARGV and die "perl $0 ipr_all6_tsv.TEprot.IPRacc > ipr_all6_tsv.TEprot.IPRacc.line\n";

my %h;
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $h{$ta[0]}{'ipr'}{$ta[11]} = $ta[12];
}
for my $id1 (sort keys %h) {
  my @k1 = sort keys %{$h{$id1}{'ipr'}};
  my @v1 = @{$h{$id1}{'ipr'}}{@k1};
  print join("\t", $id1, join(";;", @k1), join(";;", @v1))."\n";
}
