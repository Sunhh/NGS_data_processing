#!/usr/bin/perl
use strict; 
use warnings; 

!@ARGV and die "perl $0 full_table.csv > full_table.cnt\n"; 

my @out_k = qw(Complete_any Complete_1 Complete_2 Complete_high Fragmented Missing); 
print join("\t", "File", "BUSCO_gene_number", @out_k, map { "$_\%" } @out_k)."\n"; 
my @kk = qw/Complete Duplicated Fragmented Missing/; 
my %h; 
for my $fn (@ARGV) {
  open F,'<',"$fn" or die ; 
  my %h2; 
  while (<F>) {
    m!^#! and next; 
    chomp; 
    my @ta=split(/\t/, $_); 
    $h2{$ta[1]}{$ta[0]} ++; 
  }
  close F; 
  for (@kk) {
    $h2{$_} //= {}; 
  }
  $h{$fn}{'Complete_1'} = scalar(keys %{$h2{'Complete'}}); 
  $h{$fn}{'Complete_any'} = scalar(keys %{$h2{'Complete'}}) + scalar(keys %{$h2{'Duplicated'}}); 
  $h{$fn}{'Missing'} = scalar(keys %{$h2{'Missing'}}); 
  $h{$fn}{'Fragmented'} = scalar(keys %{$h2{'Fragmented'}}); 
  for my $k1 (keys %{$h2{'Duplicated'}}) {
    if ($h2{'Duplicated'}{$k1} == 2) {
      $h{$fn}{'Complete_2'} ++; 
    } else {
      $h{$fn}{'Complete_high'} ++; 
    }
  }
  $h{$fn}{'Complete_2'} //= 0; 
  $h{$fn}{'Complete_high'} //= 0; 
  my %perc; 
  my $sum = 0;
  for (qw(Complete_any Fragmented Missing)) { $sum += $h{$fn}{$_}; }
  for (@out_k) { $perc{$_} = sprintf("%0.2f", 100 * $h{$fn}{$_}/$sum); }

  print join("\t", $fn, $sum, @{$h{$fn}}{@out_k}, @perc{@out_k})."\n"; 
}

