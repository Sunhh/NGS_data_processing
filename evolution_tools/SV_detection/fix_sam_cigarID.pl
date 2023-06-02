#!/usr/bin/perl
# 5/4/2023 Query must be aligned.
# 5/30/2023 Remove ending insertions/deletions which may cause problem in NucDiff. Editing.
# 6/2/2023: Fix NM:i:\d+ when cigar is changed.

use strict;
use warnings;
use fileSunhh;
use SeqAlnSunhh;

-t and !@ARGV and die "perl $0 in.sam > restored.sam\n";

my $rev_flag = &SeqAlnSunhh::mk_flag('keep'=>'2=0,4=1');
my $fwd_flag = &SeqAlnSunhh::mk_flag('keep'=>'2=0,4=0');
my %flag2str;
for (keys %$rev_flag) {$flag2str{$_}='r';}
for (keys %$fwd_flag) {$flag2str{$_}='f';}

while (<>) {
  chomp;
  if (m!^\@!) {
    print STDOUT "$_\n";
    next;
  }
  my @ta=split(/\t/, $_);
  $ta[2] eq '*' and next;
  defined $flag2str{$ta[1]} or die "bad 4: flag: $ta[1]\n";
  if ($flag2str{$ta[1]} eq 'f') {
    # Forward alignment.
    # Read left;
    ### Fix I.
    if ($ta[5] =~ s!^(\d+)H(\d+)I!!) {
      my ($vH,$vV)=($1,$2);
      $ta[5] = ($vH+$vV)."H$ta[5]"; $ta[9] eq '*' or substr($ta[9], 0, $vV)='';
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    } elsif ($ta[5] =~ s!^(\d+)I!!) {
      my $vV = $1;
      $ta[5] = $vV."H$ta[5]"; $ta[9] eq '*' or substr($ta[9], 0, $vV)='';
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    }
    ### Fix D.
    if ($ta[5] =~ s!^(\d+)H(\d+)D!!) {
      my ($vH,$vV)=($1,$2);
      $ta[5] = $vH."H$ta[5]"; $ta[3] += $vV;
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    } elsif ($ta[5] =~ s!^(\d+)D!!) {
      my $vV = $1;
      $ta[3] += $vV;
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    }
    # Read right.
    ### Fix I.
    if ($ta[5] =~ s!(\d+)I(\d+)H$!!) {
      my ($vV,$vH)=($1,$2);
      $ta[5] = $ta[5] . ($vV+$vH) . "H"; $ta[9] eq '*' or substr($ta[9], -$vV)='';
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    } elsif ($ta[5] =~ s!(\d+)I$!!) {
      my $vV = $1;
      $ta[5] = $ta[5] . $vV . "H"; $ta[9] eq '*' or substr($ta[9], -$vV)='';
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    }
    ### Fix D.
    if ($ta[5] =~ s!(\d+)D(\d+)H$!!) {
      my ($vV,$vH)=($1,$2);
      $ta[5] = $ta[5].$vH."H";
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    } elsif ($ta[5] =~ s!(\d+)D$!!) {
      my $vV = $1;
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    }
  } else {
    # Reverse alignment.
    # Read left.
    ### Fix I.
    if ($ta[5] =~ s!^(\d+)H(\d+)I!!) {
      my ($vH,$vV)=($1,$2);
      $ta[5] = ($vH+$vV)."H$ta[5]"; $ta[9] eq '*' or substr($ta[9], 0, $vV)='';
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    } elsif ($ta[5] =~ s!^(\d+)I!!) {
      my $vV = $1;
      $ta[5] = $vV."H$ta[5]"; $ta[9] eq '*' or substr($ta[9], 0, $vV)='';
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    }
    ### Fix D.
    if ($ta[5] =~ s!^(\d+)H(\d+)D!!) {
      my ($vH,$vV)=($1,$2);
      $ta[5] = $vH."H$ta[5]"; $ta[3] += $vV;
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    } elsif ($ta[5] =~ s!^(\d+)D!!) {
      my $vV = $1;
      $ta[3] += $vV;
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    }
    # Read right.
    ### Fix I.
    if ($ta[5] =~ s!(\d+)I(\d+)H$!!) {
      my ($vV,$vH)=($1,$2);
      $ta[5] = ($vV+$vH)."H$ta[5]"; $ta[9] eq '*' or substr($ta[9], -$vV)='';
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    } elsif ($ta[5] =~ s!(\d+)I$!!) {
      my $vV = $1;
      $ta[5] = $vV."H$ta[5]"; $ta[9] eq '*' or substr($ta[9], -$vV)='';
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    }
    ### Fix D.
    if ($ta[5] =~ s!(\d+)D(\d+)H$!!) {
      my ($vV,$vH)=($1,$2);
      $ta[5] = $ta[5] . $vH . "H";
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    } elsif ($ta[5] =~ s!(\d+)D$!!) {
      my $vV = $1;
      for my $tb (@ta) { $tb =~ m!^NM:i:(\d+)$! and $tb = "NM:i:".($1-$vV); }
    }
  }
  print join("\t",@ta)."\n";
}

