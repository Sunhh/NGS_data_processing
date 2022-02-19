#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0 in.fa.N50 in_2.fa.N50 > out.table\n";

print join("\t", qw/Accession asm_size asm_num asm_N50 asm_N90/)."\n";
my @oKey = qw/ttl_bp ttl_num N50 N90/;
for my $fn (@ARGV) {
  my %h;
  open F,'<',"$fn" or die;
  while (<F>) {
    m!^Total sequences number.*:\s*(\d+)! and $h{'ttl_num'} = $1;
    m!^Total sequences bp \(ATGC\):\s*(\d+)! and $h{'ttl_bp'} = $1;
    m!^N50.+:\s*(\d+)! and $h{'N50'} = $1;
    m!^N90.+:\s*(\d+)! and $h{'N90'} = $1;
  }
  close F;
  print join("\t", $fn, @h{@oKey})."\n";
}

