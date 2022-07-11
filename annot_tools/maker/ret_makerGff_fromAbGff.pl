#!/usr/bin/perl
# [4/20/2022] Transform the maker ab-init gff3 gene models (composed of match/match_part elements) to normal gene models (gene/mRNA/exon/CDS/...)
#   Need to set up a proper environment for maker software.
use strict;
use warnings;
use LogInforSunhh;

!@ARGV and die "perl $0 out_prefix genome.fa in_maker_ab.gff3\n";

my $opref = shift;
my $genomeFa = shift;
my $abGff = shift;

my $idxPref = $genomeFa;
$idxPref =~ s!^\S+/!!;
$idxPref =~ s!\.fa(sta)?$!!i;

my $exe_maker = "/data/Sunhh/src/Support/mpich/maker_mpich/install/bin/mpiexec -n 20 /data/Sunhh/src/Annotation/maker/maker.3.01.03/bin/maker ";
my $pl_rmMkFa = "perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/rm_maker_fasta.pl ";

&runCmd("maker -CTL");
open F1,'<',"maker_opts.ctl" or die "maker_opts.ctl\n";
open O1,'>',"maker_opts.ab.ctl" or die "maker_opts.ab.ctl\n";
while (<F1>) {
  s!^genome=\s!genome=$genomeFa !
    or s!^model_org=\S+!model_org= !
    or s!^repeat_protein=\S+!repeat_protein= !
    or s!^pred_gff=\s!pred_gff=$abGff !
    or s!^keep_preds=0!keep_preds=1!;
  print O1 $_;
}
close O1;
close F1;

&runCmd("$exe_maker maker_opts.ab.ctl maker_bopts.ctl maker_exe.ctl 1> maker.ab-log 2> maker.ab-err");
&runCmd("gff3_merge -d $idxPref.maker.output/${idxPref}_master_datastore_index.log -g -o ${opref}_maker.gff3");
&runCmd("$pl_rmMkFa ${opref}_maker.gff3 | grep -v ^# > ${opref}_maker_noFa.gff3");
&runCmd("fasta_merge -d $idxPref.maker.output/${idxPref}_master_datastore_index.log -o $opref");



