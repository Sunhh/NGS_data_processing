#!/usr/bin/perl
use strict;
use warnings;


while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  if ($ta[0] eq 'Tags') {
    $ta[10] eq 'TestSet Sequences' or die "Bad columns\n";
    print join("\t", @ta[1..5], 'DEGs annotated', 'Genes annotated', 'Total DEGs', 'Total genes', "Genes")."\n";
    next;
  }
  my $ttl_geneN = $ta[6]+$ta[7]+$ta[8]+$ta[9];
  my $ttl_degN  = $ta[6]+$ta[8];
  print join("\t", @ta[1..5], $ta[6], $ta[7], $ttl_degN, $ttl_geneN, $ta[10])."\n";
}

# 0   Tags
# 1   GO ID
# 2   GO Name : description of the category;
# 3   GO Category: category of the function;
# 4   FDR
# 5   P-Value
# 6   Nr Test : number of transcripts in the sample with that function
# 7   Nr Reference: number of transcripts in the reference transcriptome with that function
# 8   Non Annot Test: number of transcripts without that function in the sample.
# 9   Non Annot Reference: number of transcripts without that function in the reference.
# 10  TestSet Sequences

