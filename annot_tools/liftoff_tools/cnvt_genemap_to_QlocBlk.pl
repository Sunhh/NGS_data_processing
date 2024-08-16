#!/usr/bin/perl
# [3/28/2022] Retrieve CDS sequences if they have not been predicted yet.
# 10/4/2023: Allow gene-location IDs without group tag.
# 10/7/2023: Output table instead of FASTA file.
use strict;
use warnings;
use fastaSunhh;
my $fs_obj = fastaSunhh->new();

!@ARGV and die "perl $0 genomic.fa  output/mapCDS.CLpan.to.CMpan.tbl.gene_map > output/mapCDS.CLpan.to.CMpan.tbl.gene_map.Qcds.blk\n";

my $fnFas = shift;
my %seq = %{ $fs_obj->save_seq_to_hash( 'faFile' => $fnFas ) };
for (keys %seq) { $seq{$_}{'seq'} =~ s!\s!!g; $seq{$_}{'len'} = length($seq{$_}{'seq'}); }

my %h;
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'ori_R_ID' and next;
  $ta[1] =~ m!^(?:[^\s:]+:)?\S+:\d+\-\d+:[+-]$! or next;
  defined $h{$ta[1]} and next;
  $h{$ta[1]} = 1;
  my $chr_id  = $ta[4];
  my $chr_str = $ta[7];
  my $cds_seq = '';
  defined $seq{$chr_id} or die "$chr_id\n";
  my @segP = split(/;/, $ta[17]);
  $chr_str eq '-' and @segP = reverse(@segP);
  for my $tb (@segP) {
    $tb =~ m!^(\d+)\,(\d+)$! or die "$tb\n";
    $cds_seq .= substr($seq{$chr_id}{'seq'}, $1-1, $2-$1+1);
  }
  $chr_str eq '-' and &fastaSunhh::rcSeq(\$cds_seq, 'rc');
  print STDOUT join("\t", $ta[1], $chr_id, $chr_str, $ta[17], $cds_seq, length($cds_seq))."\n";
}

