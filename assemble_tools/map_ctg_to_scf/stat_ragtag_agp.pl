#!/usr/bin/perl
use strict;
use warnings;
use LogInforSunhh;
use mathSunhh;

my $htxt = <<HH;
################################################################################
perl $0 '^Cla97Chr' BC19/noRedNoContRef97/ragtag.scaffold.agp > BC19/noRedNoContRef97/ragtag.scaffold.agp.stats

 '^Cla97Chr' is a perl regular expression pattern to define a chromosome ID.


HH

!@ARGV and &LogInforSunhh::usage($htxt);

my $ptn = shift; $ptn = qr/$ptn/s;

my @okeys = qw/all_seqN all_bp placed_seqN placed_bp unplaced_seqN unplaced_bp gap_bp gap_seqN/;
print STDOUT join("\t", qw/filename/, @okeys)."\n";
for my $fn (@ARGV) {
  open F,'<',"$fn" or die;
  my %h;
  while (<F>) {
    m!^\s*($|#)! and next;
    chomp;
    my @ta=split(/\t/, $_);
    if ($ta[4] eq 'U' or $ta[4] eq 'N') {
      $h{'gap_seqN'} ++;
      $h{'gap_bp'} += $ta[5];
      next;
    } elsif ($ta[4] ne 'W') {
      &stopErr("[Err] Unknown ta[4] [$ta[4]] in line: $_\n");
    }
    $h{'all_seqID'}{$ta[5]} ++;
    $h{'all_bp'} += $ta[7]-$ta[6]+1;
    if ($ta[0] =~ m!$ptn!o) {
      $h{'placed_seqID'}{$ta[5]} ++;
      $h{'placed_bp'} += $ta[7]-$ta[6]+1;
    } else {
      $h{'unplaced_seqID'}{$ta[5]} ++;
      $h{'unplaced_bp'} += $ta[7]-$ta[6]+1;
    }
  }
  close F;
  $h{'all_seqN'} = scalar(keys %{$h{'all_seqID'}});
  $h{'placed_seqN'} = scalar(keys %{$h{'placed_seqID'}});
  $h{'unplaced_seqN'} = scalar(keys %{$h{'unplaced_seqID'}});
  print STDOUT join("\t", $fn, @h{@okeys})."\n";
}

