#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0 CL comb.grp2.novl_loc.tbl output/CL.mat.pav_rd1_cov50.good4pav > output/comb.grp2.CL_PAV\n";

my $tag = shift;
my $fn_gene2grp = shift;
my $fn_mat = shift;

open F1,'<',"$fn_gene2grp" or die;
my (%gene2grp, %grpInfo);
while (<F1>) {
  chomp;
  my @ta=split(/\t/, $_);
  $gene2grp{$ta[0]} = $ta[1];
  $grpInfo{$ta[1]}{'size'} ++;
}
close F1;

my $accN;
open F2,'<',"$fn_mat" or die;
my @h1;
while (<F2>) {
  chomp;
  my @ta=split(/\t/, $_);
  if ($ta[0] eq "geneID") {
    @h1 = @ta;
    $h1[0] = "groupID";
    next;
  }
  defined $gene2grp{$ta[0]} or next;
  $accN //= $#ta;
  $accN < $#ta and die "[Err] The line length is different in file [$fn_mat]!\n";
  my $grpID = $gene2grp{$ta[0]};
  for (my $i=1; $i<@ta; $i++) {
    $grpInfo{$grpID}{'pav'}[$i-1] //= 0;
    $ta[$i] eq "P" or next;
    $grpInfo{$grpID}{'pav'}[$i-1] ++;
  }
}
close F2;
print STDOUT join("\t", @h1)."\n";
for my $grpID (sort keys %grpInfo) {
  for (my $i=0; $i<=$accN; $i++) {
    $grpInfo{$grpID}{'pav'}[$i] //= 0;
  }
  print STDOUT join("\t", $grpID, map { ($_ >= 1) ? "P" : "A" ; } @{$grpInfo{$grpID}{'pav'}})."\n";
}

