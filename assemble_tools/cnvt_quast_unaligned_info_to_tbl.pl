#!/usr/bin/perl
use strict; 
use warnings;

-t and !@ARGV and die "perl $0 ref_hc0/contigs_reports/contigs_report_hf2-noRed.unaligned.info | less -S\n"; 


print join("\t", qw/Contig Total_length Unaligned_length Unaligned_type unAln_start unAln_end segment_len/)."\n"; 
while (<>) {
  chomp;
  my @ta=split(/\t/, $_); 
  $ta[1] eq 'Total_length' and next; 
  my @tb = split(/,/, $ta[4]); 
  for my $tc (@tb) {
    $tc =~ m!^(\d+)\-(\d+)$! or die "$tc in $_\n"; 
    my ($s,$e) = ($1, $2); 
    print join("\t", @ta[0,1,2,3], $s, $e, $e-$s+1)."\n"; 
  }
}

