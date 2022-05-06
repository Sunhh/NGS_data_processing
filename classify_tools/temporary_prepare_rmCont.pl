#!/usr/bin/perl
# [4/9/2022] Need to convert this as pipeline some day later for repeat works!!!
# [5/5/2022] Fit diamond blastx to Nr database.
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "bn_task:s",
  "bn_db:s", # nt
  "help!"
);


my $htxt = <<HH;
###################################################################
perl $0  max_length  in_ctg.fa  out_prefix

-bn_task   [megablast] Could be blastx for Nr alignment.
-bn_db     [nt]

* Please set up the correct BLASTDB environment for Nt database.
* Please set up an Nt database together with taxonomy information.
* Please allow correct version of blastn program in PATH.
* Some other scripts in NGS_data_processing are required.

HH

$opts{'bn_task'} //= 'megablast';
if ($opts{'bn_task'} eq 'blastx') {
  $opts{'bn_db'} //= '/data/share/database/db_blast/NCBI/diamond_db/nr_20211207_wtax'; # This is on panda server.
} else {
  $opts{'bn_db'} //= 'nt';
}

!@ARGV and &LogInforSunhh::usage($htxt);

my $max_len = shift;
my $ctgFa   = shift;
my $opref   = shift;

&runCmd("deal_fasta.pl -keep_len 0-$max_len $ctgFa > $opref.tochk.fa");
if ($opts{'bn_task'} eq 'blastx') {
  &runCmd("perl /home/Sunhh/tools/github/NGS_data_processing/classify_tools/run_seg_blastn.pl  $opref.tochk2Nt  1000   $opref.tochk.fa  $opts{'bn_db'}  -bn_task $opts{'bn_task'}");
} else {
  &runCmd("perl /home/Sunhh/tools/github/NGS_data_processing/classify_tools/run_seg_blastn.pl  $opref.tochk2Nt  10000  $opref.tochk.fa  $opts{'bn_db'}  -bn_task $opts{'bn_task'}");
}
{
  # I'd like to keep rDNA as good assembly for contamination removal. 
  # However, it might be worth another try to treat them as contamination for chromosome anchoring.
  my $cmd = "perl /home/Sunhh/tools/github/NGS_data_processing/classify_tools/classify_region_byBn6.pl ";
  $cmd .= " $opref.tochk2Nt.bn6  -joinInEx $opref.tochk2Nt.bn6.jnInEx ";
  $cmd .= " -InList Eukaryota:Satellite:rDNA:Chloroplast:Mitochondrion:Plastid ";
  $cmd .= " -ExList NA:Viruses:Bacteria:Archaea ";
  $cmd .= " > $opref.tochk2Nt.bn6.class ";
  &runCmd($cmd);
}
&runCmd("perl /home/Sunhh/tools/github/NGS_data_processing/classify_tools/get_Ex_region.pl $opref.tochk2Nt.bn6.class 1> $opref.tochk2Nt.bn6.class.sep 2> $opref.tochk2Nt.bn6.class.jn");
if ($opts{'bn_task'} eq 'blastx') {
  &runCmd("perl /home/Sunhh/tools/github/NGS_data_processing/classify_tools/recog_organelle_rDNA_from_classJn.pl $opref.tochk2Nt.bn6.class.jn -for_bx6 > $opref.tochk2Nt.bn6.class.jn.s1");
} else {
  &runCmd("perl /home/Sunhh/tools/github/NGS_data_processing/classify_tools/recog_organelle_rDNA_from_classJn.pl $opref.tochk2Nt.bn6.class.jn > $opref.tochk2Nt.bn6.class.jn.s1");
}

# In order to detect rDNA and organelle genomes.
{
  my $cmd = "perl /home/Sunhh/tools/github/NGS_data_processing/classify_tools/classify_region_byBn6.pl ";
  $cmd .= " $opref.tochk2Nt.bn6  -joinInEx $opref.tochk2Nt.bn6.highCopy.jnInEx ";
  $cmd .= " -InList Eukaryota:NA:Viruses:Bacteria:Archaea ";
  $cmd .= " -ExList Chloroplast:Mitochondrion:Plastid:Satellite:rDNA ";
  $cmd .= " > $opref.tochk2Nt.bn6.highCopy.class ";
  &runCmd($cmd);
}
&runCmd("perl /home/Sunhh/tools/github/NGS_data_processing/classify_tools/get_Ex_region.pl $opref.tochk2Nt.bn6.highCopy.class 1> $opref.tochk2Nt.bn6.highCopy.class.sep 2> $opref.tochk2Nt.bn6.highCopy.class.jn");
if ($opts{'bn_task'} eq 'blastx') {
  &runCmd("perl /home/Sunhh/tools/github/NGS_data_processing/classify_tools/recog_organelle_rDNA_from_classJn.pl $opref.tochk2Nt.bn6.highCopy.class.jn -for_bx6 > $opref.tochk2Nt.bn6.highCopy.class.jn.s1");
} else {
  &runCmd("perl /home/Sunhh/tools/github/NGS_data_processing/classify_tools/recog_organelle_rDNA_from_classJn.pl $opref.tochk2Nt.bn6.highCopy.class.jn > $opref.tochk2Nt.bn6.highCopy.class.jn.s1");
}

# deal_fasta.pl -keep_len 0-1000000 ../db/C31.hf2a.noRed.fa > C31.tochk.fa
# blastn -query C31.tochk.fa -out C31.tochk.fa.toNt.bn6 \
#  -db nt -evalue 1e-5 -num_threads 50 -max_hsps 50 -max_target_seqs 50 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand staxids sscinames sskingdoms stitle'
# perl /home/Sunhh/tools/github/NGS_data_processing/classify_tools/classify_region_byBn6.pl \
#  C31.tochk.fa.toNt.bn6 \
#  -joinInEx   C31.tochk.fa.toNt.bn6.jnInEx \
#  -InList     Eukaryota:Satellite \
#  -ExList     NA:Viruses:Bacteria:Archaea:rDNA:Chloroplast:Mitochondrion:Plastid \
# > C31.tochk.fa.toNt.bn6.class
# perl /home/Sunhh/tools/github/NGS_data_processing/classify_tools/get_Ex_region.pl C31.tochk.fa.toNt.bn6.class 1> C31.tochk.fa.toNt.bn6.class.sep 2> C31.tochk.fa.toNt.bn6.class.jn
# awk 'NR != 1 { print $5/$2*100"\t"$0 }' C31.tochk.fa.toNt.bn6.class.jn > C31.tochk.fa.toNt.bn6.class.jn.1
# cat C31.tochk.fa.toNt.bn6.class.jn.1 | grep rDNA    | awk '{print $1"\t"$2"\trDNA"}'      >  C31.tochk.fa.toNt.bn6.class.jn.2
# cat C31.tochk.fa.toNt.bn6.class.jn.1 | grep -v rDNA | awk '{print $1"\t"$2"\torganelle"}' >> C31.tochk.fa.toNt.bn6.class.jn.2


