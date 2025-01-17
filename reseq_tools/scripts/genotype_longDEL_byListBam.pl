#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;
use LogInforSunhh;

!@ARGV and die "perl $0 list.long_deletions list.sample_bam out_prefix\n";

my $inFnDelList = shift;
my $inFnSampleBam = shift;
my $opref = shift;

my $bin_path = &fileSunhh::_abs_path($0);
$bin_path = &fileSunhh::_dirname($bin_path);
my $pl_cntRd = "perl $bin_path/cntRd_spanJunctionSite_inBam.pl";
my $pl_runCmd = "run_cmd_in_batch.pl";
my $cpuN = 60;

my @del_list = &fileSunhh::load_tabFile($inFnDelList);
my @acc2bam = map { [$_->[0], $_->[1]] } &fileSunhh::load_tabFile($inFnSampleBam);
open O1,'>',"$opref-melt.tab" or die "$opref.tab $!\n";
for my $d1 (@del_list) {
  my ($del_id1, $del_chr, $del_s, $del_e) = @$d1;
  my ($s1,$e1) = ($del_s-11, $del_s+10); $s1 < 1 and $s1 = 1; # left boundary.
  my ($s2,$e2) = ($del_e-10, $del_e+11); $s2 < 1 and $s2 = 1; # right boundary.
  my $wdir = &fileSunhh::new_tmp_dir('create'=>1);
  for my $x1 (@acc2bam) {
    my ($acc1, $bamFn) = @$x1;
    &fileSunhh::write2file("$wdir/cmd1_$del_id1", "$pl_cntRd -in_bam $bamFn -target_loc $del_chr:$s1-$e1 -positive_control_loc $del_chr:$s2-$e2 > $wdir/out.$acc1\n", ">>");
  }
  &runCmd("$pl_runCmd -cpuN $cpuN $wdir/cmd1_$del_id1 > $wdir/cmd1_$del_id1.log");
  for my $x1 (@acc2bam) {
    my ($acc1, $bamFn) = @$x1;
    open F1,'<',"$wdir/out.$acc1" or die;
    while (<F1>) {
      m!^\s*$! and next; chomp; my ($fn1, $cnt1, $cnt2)=split(/\t/, $_);
      my $t_geno = "1/1"; $cnt1 > 0 and $cnt2 > 0 and $t_geno = "0/0";
      print O1 join("\t", $del_id1, $acc1, $t_geno, $cnt1, $cnt2)."\n";
    }
    close F1;
  }
  &fileSunhh::_rmtree($wdir);
}
close O1;

