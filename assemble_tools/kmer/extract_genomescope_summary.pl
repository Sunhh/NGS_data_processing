#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;

!@ARGV and die "perl $0 out.of.genomescope/BC19/summary.txt out.of.genomescope/BC20/summary.txt > summary_gscope.tbl\n";

my @outKey = qw(filename k k_depth max_genome_hap max_hetePerc max_genome_uni max_genome_rep);
print STDOUT join("\t", qw/filename kmer kmer_depth genome_size hete% unique_region repeat_region/)."\n";
for my $fn (@ARGV) {
  open F,'<',"$fn" or die;
  my %h;
  $h{'filename'} = $fn;
  while (<F>) {
    chomp;
    if (m!^k\s*=\s*(\d+)\s*$!) {
      $h{'k'} = $1;
    } elsif (m!^Genome Haploid Length\s+[\d,]+\s*bp\s+([\d,]+)\s*bp!) {
      $h{'max_genome_hap'} = $1;
      $h{'max_genome_hap'} =~ s!,!!g;
    } elsif (m!^Genome Repeat Length\s+[\d,]+\s*bp\s+([\d,]+)\s*bp!) {
      $h{'max_genome_rep'} = $1;
      $h{'max_genome_rep'} =~ s!,!!g;
    } elsif (m!^Genome Unique Length\s+[\d,]+\s*bp\s+([\d,]+)\s*bp!) {
      $h{'max_genome_uni'} = $1;
      $h{'max_genome_uni'} =~ s!,!!g;
    } elsif (m!^Heterozygosity\s+[\d.]+\%\s+([\d.]+)\%!) {
      $h{'max_hetePerc'} = $1;
    }
  }
  close F;
  my $dirName = &fileSunhh::_dirname( &fileSunhh::_abs_path($fn) );
  if (-e "$dirName/model.txt") {
    open F2,'<',"$dirName/model.txt" or die "[Err] Failed to open file $dirName/model.txt\n";
    while (<F2>) {
      chomp;
      m!^kmercov\s+(\S+)! or next;
      $h{'k_depth'} = $1 * 2;
    }
    close F2;
  }
  for (@outKey) {
    $h{$_} //= "N/A";
  }
  print STDOUT join("\t", @h{@outKey})."\n";
}

# Sunhh@swift:/data/Sunhh/wmhifi/genome_size$ less -S out.of.genomescope/BC19/summary.txt 
# GenomeScope version 1.0
# k = 81
# 
# property                      min               max               
# Heterozygosity                0.0523306%        0.0526906%        
# Genome Haploid Length         414,606,094 bp    414,694,212 bp    
# Genome Repeat Length          79,935,019 bp     79,952,008 bp     
# Genome Unique Length          334,671,075 bp    334,742,204 bp    
# Model Fit                     97.7185%          98.7888%          
# Read Error Rate               0.226308%         0.226308%         

