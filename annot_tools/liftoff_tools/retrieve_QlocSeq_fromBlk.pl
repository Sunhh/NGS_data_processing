#!/usr/bin/perl
# [3/28/2022] Retrieve CDS sequences if they have not been predicted yet.
use strict;
use warnings;
use fastaSunhh;
my $fs_obj = fastaSunhh->new();

!@ARGV and die "perl $0 genomic.fa  output/Qcds.CX.blk.CL > output/Qcds.CX.blk.CL.fa\n";

my $fnFas = shift;
my %seq = %{ $fs_obj->save_seq_to_hash( 'faFile' => $fnFas ) };
for (keys %seq) { $seq{$_}{'seq'} =~ s!\s!!g; $seq{$_}{'len'} = length($seq{$_}{'seq'}); }

my %h;
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  my $ele_id = $ta[0];
  $ele_id eq 'Q_ID' and next;
  # $ta[0] =~ m!^C\S:\S+:\d+\-\d+:[+-]$! or die "$_\n";
  defined $h{$ta[0]} and next;
  $h{$ele_id} = 1;
  my $chr_id  = $ta[1];
  my $chr_str = $ta[2];
  my $blks    = $ta[3];
  my $cds_seq = '';
  defined $seq{$chr_id} or die "$chr_id\n";
  for my $tb (split(/;/, $blks)) {
    $tb =~ m!^(\d+)\,(\d+)$! or die "$tb\n";
    $cds_seq .= substr($seq{$chr_id}{'seq'}, $1-1, $2-$1+1);
  }
  if ($chr_str eq '-') {
    &fastaSunhh::rcSeq(\$cds_seq, 'rc');
  }
  print STDOUT ">$ele_id\n$cds_seq\n";
}

