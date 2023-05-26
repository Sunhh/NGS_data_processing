#!/usr/bin/perl
# 5/4/2023 Query must be aligned.
use strict;
use warnings;
use fileSunhh;
use SeqAlnSunhh;

!@ARGV and die "perl $0 q.fa.kl r.fa.kl in.sam > restored.sam\n";

my $klQ = shift;
my $klR = shift;

my $rev_flag = &SeqAlnSunhh::mk_flag('keep'=>'2=0,4=1');
my $fwd_flag = &SeqAlnSunhh::mk_flag('keep'=>'2=0,4=0');
my %flag2str;
for (keys %$rev_flag) {$flag2str{$_}='r';}
for (keys %$fwd_flag) {$flag2str{$_}='f';}

my %qKL;
for my $l0 (&fileSunhh::load_tabFile($klQ)) {
  $l0->[0] eq 'key' and next;
  $qKL{'id2len'}{$l0->[0]} = $l0->[1];
  push(@{$qKL{'kl'}}, @{$l0}[0,1]);
}
my %rKL;
for my $l0 (&fileSunhh::load_tabFile($klR)) {
  $l0->[0] eq 'key' and next;
  $rKL{'id2len'}{$l0->[0]} = $l0->[1];
  push(@{$rKL{'kl'}}, @{$l0}[0,1]);
  print STDOUT join("\t", '@SQ', "SN:$l0->[0]", "LN:$l0->[1]")."\n";
}

my %rsub2whole;
my %qsub2whole;
while (<>) {
  chomp;
  if (m!^\@SQ\t!) {
    m!^\@SQ\tSN:((\S+)_(\d+)_(\d+))\tLN:\d+$! or die "bad 1: $_\n";
    my ($subID, $wholeID, $subS, $subE) = ($1, $2, $3, $4);
    defined $rKL{'id2len'}{$wholeID} or die "bad 2: $_\n";
    $rsub2whole{$subID} = [$wholeID, $subS, $subE];
    next;
  } elsif (m!^\@PG\t!) {
    next;
  } elsif (m!^\@!) {
    print STDOUT "$_\n";
    next;
  }
  my @ta=split(/\t/, $_);
  $ta[2] eq '*' and next;
  defined $rsub2whole{$ta[2]} or die "bad 3: RID: $ta[2] in line: $_\n";
  defined $flag2str{$ta[1]} or die "bad 4: flag: $ta[1]\n";
  unless (defined $qsub2whole{$ta[0]}) {
    $ta[0] =~ m!^(\S+)_(\d+)_(\d+)$! or die "bad 5: QID: $ta[0]\n";
    $qsub2whole{$ta[0]} = [$1, $2, $3];
  }
  my ($qWID, $qWS, $qWE) = @{$qsub2whole{$ta[0]}};
  my ($rWID, $rWS, $rWE) = @{$rsub2whole{$ta[2]}};
  if ($flag2str{$ta[1]} eq 'f') {
    $ta[3] = $ta[3] + $rWS - 1;
    if ($ta[5] =~ s!^(\d+)S!!) {
      $ta[5] = ($1+$qWS-1)."H$ta[5]";
      $ta[9] eq '*' or substr($ta[9], 0, $1) = '';
    } elsif ($ta[5] =~ s!^(\d+)H!!) {
      $ta[5] = ($1+$qWS-1)."H$ta[5]";
    } elsif ($qWS-1 > 0) {
      $ta[5] = ($qWS-1)."H$ta[5]";
    }
    if ($ta[5] =~ s!(\d+)S$!!) {
      $ta[5] = $ta[5].($1+$qKL{'id2len'}{$qWID}-$qWE)."H";
      $ta[9] eq '*' or substr($ta[9], -$1)='';
    } elsif ($ta[5] =~ s!(\d+)H$!!) {
      $ta[5] = $ta[5].($1+$qKL{'id2len'}{$qWID}-$qWE)."H";
    } elsif ($qKL{'id2len'}{$qWID}-$qWE > 0) {
      $ta[5] = $ta[5].($qKL{'id2len'}{$qWID}-$qWE)."H";
    }
  } else {
    $ta[3] = $ta[3] + $rWS - 1;
    if ($ta[5] =~ s!^(\d+)S!!) {
      $ta[5] = ($1+$qKL{'id2len'}{$qWID}-$qWE)."H$ta[5]";
      $ta[9] eq '*' or substr($ta[9], 0, $1) ='';
    } elsif ($ta[5] =~ s!^(\d+)H!!) {
      $ta[5] = ($1+$qKL{'id2len'}{$qWID}-$qWE)."H$ta[5]";
    } elsif ($qKL{'id2len'}{$qWID}-$qWE > 0) {
      $ta[5] = ($qKL{'id2len'}{$qWID}-$qWE)."H$ta[5]";
    }
    if ($ta[5] =~ s!(\d+)S$!!) {
      $ta[5] = $ta[5].($1+$qWS-1)."H";
      $ta[9] eq '*' or substr($ta[9], -$1) ='';
    } elsif ($ta[5] =~ s!(\d+)H$!!) {
      $ta[5] = $ta[5].($1+$qWS-1)."H";
    } elsif ($qWS-1 > 0) {
      $ta[5] = $ta[5].($qWS-1)."H";
    }
  }
  $ta[0] = $qWID;
  $ta[2] = $rWID;
  print join("\t",@ta)."\n";
}

