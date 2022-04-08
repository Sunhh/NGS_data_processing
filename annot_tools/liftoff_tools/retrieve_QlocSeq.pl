#!/usr/bin/perl
# [3/28/2022] Retrieve CDS sequences if they have not been predicted yet.
use strict;
use warnings;
use fastaSunhh;
my $fs_obj = fastaSunhh->new();

!@ARGV and die "perl $0 genomic.fa  output/mapCDS.CLpan.to.CMpan.tbl.gene_map > output/mapCDS.CLpan.to.CMpan.tbl.gene_map.Qcds.fa\n";

my $fnFas = shift;
my %seq = %{ $fs_obj->save_seq_to_hash( 'faFile' => $fnFas ) };
for (keys %seq) { $seq{$_}{'seq'} =~ s!\s!!g; $seq{$_}{'len'} = length($seq{$_}{'seq'}); }

my %h;
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'ori_R_ID' and next;
  $ta[1] =~ m!^C\S:\S+:\d+\-\d+:[+-]$! or next;
  defined $h{$ta[1]} and next;
  $h{$ta[1]} = 1;
  my $chr_id  = $ta[4];
  my $chr_str = $ta[7];
  my $cds_seq = '';
  defined $seq{$chr_id} or die "$chr_id\n";
  for my $tb (split(/;/, $ta[17])) {
    $tb =~ m!^(\d+)\,(\d+)$! or die "$tb\n";
    $cds_seq .= substr($seq{$chr_id}{'seq'}, $1-1, $2-$1+1);
  }
  if ($chr_str eq '-') {
    &fastaSunhh::rcSeq(\$cds_seq, 'rc');
  }
  print STDOUT ">$ta[1]\n$cds_seq\n";
}

