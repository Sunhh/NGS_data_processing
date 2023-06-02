#!/usr/bin/perl
# 5/4/2023 Query must be aligned.
# 5/30/2023 Remove ending insertions/deletions which may cause problem in NucDiff. Editing.

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
      $ta[5] = ($1+$2)."H$ta[5]";
      $ta[9] eq '*' or substr($ta[9], 0, $2)='';
    } elsif ($ta[5] =~ s!^(\d+)I!!) {
      $ta[5] = $1."H$ta[5]";
      $ta[9] eq '*' or substr($ta[9], 0, $1)='';
    }
    ### Fix D.
    if ($ta[5] =~ s!^(\d+)H(\d+)D!!) {
      $ta[5] = $1."H$ta[5]";
      $ta[3] += $2;
    } elsif ($ta[5] =~ s!^(\d+)D!!) {
      $ta[3] += $1;
    }
    # Read right.
    ### Fix I.
    if ($ta[5] =~ s!(\d+)I(\d+)H$!!) {
      $ta[5] = $ta[5] . ($1+$2) . "H";
      $ta[9] eq '*' or substr($ta[9], -$1)='';
    } elsif ($ta[5] =~ s!(\d+)I$!!) {
      $ta[5] = $ta[5] . $1 . "H";
      $ta[9] eq '*' or substr($ta[9], -$1)='';
    }
    ### Fix D.
    if ($ta[5] =~ s!(\d+)D(\d+)H$!!) {
      $ta[5] = $ta[5].$2."H";
    } elsif ($ta[5] =~ s!(\d+)D$!!) {
      ;
    }
  } else {
    # Reverse alignment.
    # Read left.
    ### Fix I.
    if ($ta[5] =~ s!^(\d+)H(\d+)I!!) {
      $ta[5] = ($1+$2)."H$ta[5]";
      $ta[9] eq '*' or substr($ta[9], 0, $2)='';
    } elsif ($ta[5] =~ s!^(\d+)I!!) {
      $ta[5] = $1."H$ta[5]";
      $ta[9] eq '*' or substr($ta[9], 0, $1)='';
    }
    ### Fix D.
    if ($ta[5] =~ s!^(\d+)H(\d+)D!!) {
      $ta[5] = $1."H$ta[5]";
      $ta[3] += $2;
    } elsif ($ta[5] =~ s!^(\d+)D!!) {
      $ta[3] += $2;
    }
    # Read right.
    ### Fix I.
    if ($ta[5] =~ s!(\d+)I(\d+)H$!!) {
      $ta[5] = ($1+$2)."H$ta[5]";
      $ta[9] eq '*' or substr($ta[9], -$1)='';
    } elsif ($ta[5] =~ s!(\d+)I$!!) {
      $ta[5] = $1."H$ta[5]";
      $ta[9] eq '*' or substr($ta[9], -$1)='';
    }
    ### Fix D.
    if ($ta[5] =~ s!(\d+)D(\d+)H$!!) {
      $ta[5] = $ta[5] . $2 . "H";
    } elsif ($ta[5] =~ s!(\d+)D$!!) {
      ;
    }
  }
  print join("\t",@ta)."\n";
}

