#!/usr/bin/perl
# 20211122: The NM number doesn't include the part in 'H' or 'S' regions. Process output of sam2pairwise ; 
use strict; 
use warnings; 

my %blk; 
$blk{'wind'} = 50e3; 
$blk{'curr'} = 0; 

print join("\t", qw/qname flag sname spos mapq sStart sEnd identity qGap sGap subs same/)."\n"; 

my @rec; 
while ($_ = <>) {
  chomp($_); $rec[0] = $_; 
  $_ = <>; chomp($_); $rec[1] = $_; 
  $_ = <>; chomp($_); $rec[2] = $_; 
  $_ = <>; chomp($_); $rec[3] = $_; 
  $rec[1] =~ m!^\*! and next; 
  $rec[1] =~ m!^\s*$! and next; 
  # $blk{'qS'} = 1; # Need to resolve +/- strand. 
  my @ta=split(/\t/, $rec[0]); 
  $blk{'sS'} = $ta[3]; 
  $blk{'sE'} = $blk{'sS'}; 
  if ($ta[5] =~ m!^(\d+)S!) {
    # $blk{'qS'} += $1; 
    $rec[1] = substr($rec[1], $1); 
    $rec[2] = substr($rec[2], $1); 
    $rec[3] = substr($rec[3], $1); 
  }
  if ($ta[5] =~ m!(\d+)S$!) {
    my $l1 = length($rec[1]); 
    substr($rec[1], $l1-$1) = ""; 
    substr($rec[2], $l1-$1) = ""; 
    substr($rec[3], $l1-$1) = ""; 
  }
  
  for (my $i=0; $i<length($rec[3]); $i++) {
    $blk{'curr'} ++; 
    my $r_1 = substr($rec[1], $i-1, 1); 
    my $r_3 = substr($rec[3], $i-1, 1); 
    if ($r_1 eq $r_3) {
      $blk{'same'} ++; 
      $blk{'sE'} ++; 
    } elsif ($r_1 eq '-') {
      $blk{'q_gap'} ++; 
      $blk{'sE'} ++; 
    } elsif ($r_3 eq '-') {
      $blk{'s_gap'} ++; 
      if ($blk{'sS'} > $blk{'sE'}) {
        $blk{'sS'} = $blk{'sE'}; 
      }
    } else {
      $blk{'diff'} ++; 
      $blk{'sE'} ++; 
    }
    if ($blk{'curr'} % $blk{'wind'} == 0) {
      for my $k1 (qw/q_gap s_gap diff/) {
        $blk{$k1} //= 0; 
      }
      my $ttl_diff = $blk{'q_gap'} + $blk{'s_gap'} + $blk{'diff'}; 
      my $ttl_ident = sprintf("%0.2f", $blk{'same'} / ($blk{'same'} + $ttl_diff) * 100 ); 
      print join("\t", @ta[0..4], $blk{'sS'}, $blk{'sE'}, $ttl_ident, $blk{'q_gap'}, $blk{'s_gap'}, $blk{'diff'}, $blk{'same'})."\n"; 
      $blk{'sS'} = $blk{'sE'}+1; 
      for my $k1 (qw/curr same q_gap s_gap diff/) {
        $blk{$k1} = 0; 
      }
    }
  }
  if ($blk{'curr'} > 0) {
    my $ttl_diff = $blk{'q_gap'} + $blk{'s_gap'} + $blk{'diff'}; 
    my $ttl_ident = sprintf("%0.2f", $blk{'same'} / ($blk{'same'} + $ttl_diff) * 100 ); 
    print join("\t", @ta[0..4], $blk{'sS'}, $blk{'sE'}, $ttl_ident, $blk{'q_gap'}, $blk{'s_gap'}, $blk{'diff'}, $blk{'same'})."\n"; 
    $blk{'sS'} = $blk{'sE'}+1; 
    for my $k1 (qw/curr same q_gap s_gap diff/) {
      $blk{$k1} = 0; 
    }
  }
}
