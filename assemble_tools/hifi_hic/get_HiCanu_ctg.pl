#!/usr/bin/perl
# Get primary contigs from HiCanu assemblies by removing bubble elements. 

use strict; 
use warnings; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 
use fileSunhh; 

-t and !@ARGV and die "perl $0 HiCanu_asm.contigs.fasta HiCanu_asm.contigs.layout.tigInfo > HiCanu_asm_ctg.fa\n"; 

my $f1 = shift; 
my $f2 = shift; 

my %faSeq = %{ $fs_obj->save_seq_to_hash('faFile'=>$f1) }; 
my @ids = sort { $faSeq{$a}{'Order'} <=> $faSeq{$b}{'Order'} } keys %faSeq; 


$ids[0] =~ m!^(tig\d+)$! or die "bad example tigID [$ids[0]]\n"; 
my $nLen = length($ids[0]) - 3;

my $ofh2 = &openFH($f2); 

while (<$ofh2>) {
  chomp;
  my @ta=split(/\t/, $_); 
  if (m!^\s*(#|$)!) {
    # print "$ta[0]\n";
    next;
  }
  $ta[3] eq 'contig' or next; 
  $ta[5] eq 'no' or next; 
  $ta[0] =~ m!^\d+$! or die "bad number [$ta[0]]\n";
  my $newID = sprintf("%0${nLen}d", $ta[0]);
  $newID = "tig$newID";
  defined $faSeq{$newID} or die "no seq for ID [$newID]: $_\n"; 
  print STDOUT ">$faSeq{$newID}{'head'}\n$faSeq{$newID}{'seq'}\n"; 
}
close($ofh2); 


