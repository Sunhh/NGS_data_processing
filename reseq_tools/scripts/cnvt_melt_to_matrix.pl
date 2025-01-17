#!/usr/bin/perl
use strict;
use warnings;

-t and !@ARGV and die "perl $0 in.melt > out.matrix\n";

# 1_608963_v75845_65251   ARO18917        1/1     14      0
# 1_608963_v75845_65251   ARO18920        1/1     45      0
# 1_608963_v75845_65251   ARO19494        1/1     24      0
# 1_608963_v75845_65251   ARO20587        1/1     16      0

my (%col1, %col2, %col_val);

while (<>) {
  my @ta=split(/\t/, $_);
  chomp(@ta);
  $col1{$ta[0]} //= $.;
  $col2{$ta[1]} //= $.;
  $col_val{$ta[0]}{$ta[1]} //= $ta[2];
}
my @arr1 = sort {$col1{$a}<=>$col1{$b}} keys %col1;
my @arr2 = sort {$col2{$a}<=>$col2{$b}} keys %col2;
print STDOUT join("\t", "Sample", @arr1)."\n";
for my $a2 (@arr2) {
  my @o = ($a2);
  for my $a1 (@arr1) { $col_val{$a1}{$a2} //= "./."; push(@o, $col_val{$a1}{$a2}); }
  print STDOUT join("\t", @o)."\n";
}

