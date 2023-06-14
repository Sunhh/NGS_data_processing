#!/usr/bin/perl
use strict;
use warnings;

my (@gene, @mrna, @ele);
my %e2r = qw(exon 1 five_prime_UTR 1.2 CDS 2 three_prime_UTR 3);
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  if ($ta[2] eq "gene") {
    scalar(@gene) > 0 and &out(@gene, @mrna, @ele);
    @gene = @mrna = @ele = ();
    push(@gene, [@ta]);
  } elsif ($ta[2] eq 'mRNA') {
    push(@mrna, [@ta]);
  } elsif ($ta[2] =~ m!^(exon|five_prime_UTR|three_prime_UTR)$!) {
    push(@ele, [@ta]);
  } elsif ($ta[2] eq 'CDS') {
    $ta[2] eq "CDS" and $ta[8] =~ s!^ID[^\s;=]+;!!;
    push(@ele, [@ta]);
  }
}
scalar(@gene) > 0 and &out(@gene, @mrna, @ele);
@gene = @mrna = @ele = ();

sub out {
  my (@a1) = @_;
  $a1[0][2] eq "gene" or die "bad 1 @{$a1[0]}\n";
  $a1[1][2] eq "mRNA" or die "bad 2 @{$a1[0]}\n";
  if ($a1[0][6] eq "+") {
    @a1[2 .. $#a1] = sort { $a->[3] <=> $b->[3] || $e2r{$a->[2]} <=> $e2r{$b->[2]} } @a1[2 .. $#a1];
  } elsif ($a1[0][6] eq "-") {
    @a1[2 .. $#a1] = sort { $b->[4] <=> $a->[4] || $e2r{$a->[2]} <=> $e2r{$b->[2]} } @a1[2 .. $#a1];
  } else {
    die "bad 3 @{$a1[0]}\n";
  }
  $a1[2][2] eq "exon" or die "bad 4a:@{$a1[2]}\nbad 4b:@{$a1[3]}\n";
  my $i = 3;
  for ($i=3;$i<@a1;$i++) {
    $a1[$i][2] eq 'five_prime_UTR' and next;
    $a1[$i][2] eq 'exon' and next;
    $a1[$i][2] eq 'CDS' or die "i=$i; @{$a1[$i]}\ni=$i-1;@{$a1[$i-1]}\n";
    last;
  }
  $a1[$i][2] eq "CDS" or die "bad 5:i=$i:@{$a1[$i]}\n";
  if ($a1[$i][7] == 1) {
    $a1[$i][7] = 2;
  } elsif ($a1[$i][7] == 2) {
    $a1[$i][7] = 1;
  }
  for my $a2 (@a1) {
    print join("\t", @$a2)."\n";
  }
}
