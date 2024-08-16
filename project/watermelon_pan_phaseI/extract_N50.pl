#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0 in.fa.N50 in_2.fa.N50 > out.table\n";

my @oKey = qw/asm_size asm_num longest asm_N25 asm_N50 asm_N90 asm_N95 asm_N99 shortest gnm_size GN90 GN50/;
print join("\t", 'accession', @oKey)."\n";
for my $fn (@ARGV) {
  my %h;
  open F,'<',"$fn" or die;
  while (<F>) {
    m!^Total sequences number.*:\s*(\d+)! and $h{'asm_num'} //= $1;
    m!^Total sequences bp \(ATGC\):\s*(\d+)! and $h{'asm_size'} //= $1;
    m!^Maximum length \(ATGC\)\s*:\s*(\d+)! and $h{'longest'} //= $1;
    m!^Minimum length \(ATGC\)\s*:\s*(\d+)! and $h{'shortest'} //= $1;
    m!^N25.+:\s*(\d+)! and $h{'asm_N25'} //= $1;
    m!^N50.+:\s*(\d+)! and $h{'asm_N50'} //= $1;
    m!^N90.+:\s*(\d+)! and $h{'asm_N90'} //= $1;
    m!^N95.+:\s*(\d+)! and $h{'asm_N95'} //= $1;
    m!^N99.+:\s*(\d+)! and $h{'asm_N99'} //= $1;
    m!^Est\. Genome size\s*:\s*(\d+)! and $h{'gnm_size'} //= $1;
    m!^NG50.+:\s*(\d+)! and $h{'GN50'} //= $1;
    m!^NG90.+:\s*(\d+)! and $h{'GN90'} //= $1;
  }
  close F;
  $h{'gnm_size'} == 0 and $h{'gnm_size'} = 'NA';
  for my $k1 (@oKey) {
    $h{$k1} //= 'NA';
  }
  print join("\t", $fn, @h{@oKey})."\n";
}

