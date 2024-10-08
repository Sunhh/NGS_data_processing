#!/usr/bin/perl
# [8/9/2022] Compress .log files.
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;

!@ARGV and die "perl $0 busco_output_dir_toRM/\n";

for my $dd (@ARGV) {
  opendir my $dir, "$dd" or die;
  my @files = grep { $_ !~ m!^\.+$! } readdir($dir);
  closedir($dir);
  for my $f1 (@files) {
    if ($f1 =~ m!^run_\S+_odb\d+$!) {
      -e "$dd/$f1/hmmer_output/" and &fileSunhh::_rmtree("$dd/$f1/hmmer_output/");
      -e "$dd/$f1/busco_sequences/" and &fileSunhh::_rmtree("$dd/$f1/busco_sequences/");
      -e "$dd/$f1/metaeuk_output/" and &fileSunhh::_rmtree("$dd/$f1/metaeuk_output/");
    }
  }
  for my $a1 (qw/busco hmmsearch_out metaeuk_out hmmsearch_err metaeuk_err/) {
    -e "$dd/logs/$a1.log" and &runCmd("bgzip -@ 10 $dd/logs/$a1.log");
  }
}

# rm -rf prot_*/run_embryophyta_odb10/hmmer_output/
# rm -rf prot_*/run_embryophyta_odb10/busco_sequences/
# rm -rf genom_*/run_embryophyta_odb10/hmmer_output/
# rm -rf genom_*/run_embryophyta_odb10/busco_sequences/
# rm -rf genom_*/run_embryophyta_odb10/metaeuk_output/
# gzip genom_*/logs/*.log



