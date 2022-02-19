#!/usr/bin/perl
use strict; 
use warnings;
use fastaSunhh;
use LogInforSunhh; 
my $fs_obj = fastaSunhh->new();


# WM1029_PI525083_NODE_16105_1584-62701370-4687   WM1029_PI525083_NODE_16105_1584-6270    100.00  3318    0       0       1       3318    1370    4687    0.0     6128    3318    4687    plus
!@ARGV and die "perl $0 source.fa  dedup.fa  > dedup2source.tab\n"; 

my $f1_source = shift; 
my $f2_dedup = shift; 

&tsmsg("[Log] Loading [$f1_source]\n");
my %seq1 = %{ $fs_obj->save_seq_to_hash( 'faFile' => $f1_source ) }; 
&tsmsg("[Log] Loading [$f2_dedup]\n");
my %seq2 = %{ $fs_obj->save_seq_to_hash( 'faFile' => $f2_dedup  ) }; 

&tsmsg("[Log] Cleaning sequences.\n");
my @seq1ID = sort keys %seq1; 
my @seq2ID = sort keys %seq2; 
for (@seq1ID) { $seq1{$_}{'seq'} =~ s!\s!!g; $seq1{$_}{'len'} = length($seq1{$_}{'seq'}); }
for (@seq2ID) { $seq2{$_}{'seq'} =~ s!\s!!g; $seq2{$_}{'len'} = length($seq2{$_}{'seq'}); }

print STDOUT join("\t", qw/dedup_ID rmcont_ID rmcont_start rmcont_end/)."\n";
my $ttl = scalar(@seq2ID); 
my $cnt2 = 0; 
for my $id2 (@seq2ID) {
  $cnt2 ++; 
  my $proc_perc = sprintf("%.2f", 100*$cnt2/$ttl); 
  &tsmsg("[Log] Processing [$id2] ($proc_perc%)\n"); 
  my $is_have = 0;
  for my $id1 (@seq1ID) {
    if ($id1 eq $id2 and $seq1{$id1}{'seq'} eq $seq2{$id2}{'seq'}) {
      print STDOUT join("\t", $id2, $id1, 1, $seq1{$id1}{'len'})."\n"; 
      $is_have = 1; 
      last; 
    } elsif ($id2 =~ m!${id1}(\d+)\-(\d+)$!) {
      my ($s, $e) = ($1, $2); 
      if ( substr($seq1{$id1}{'seq'}, $s-1, $e-$s+1) eq $seq2{$id2}{'seq'} ) {
        print STDOUT join("\t", $id2, $id1, $s, $e)."\n"; 
        $is_have = 1;
        last; 
      }
    }
  }
  if ($is_have == 0) {
    &tsmsg("[Wrn] Failed to map dedup_ID [$id2]\n"); 
  }
}


