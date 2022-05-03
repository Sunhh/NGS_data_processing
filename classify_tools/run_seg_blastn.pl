#!/usr/bin/perl
# [4/20/2022] Because the input contigs are super long, it may lower the sensitivity by setting lower -max_hsp and -max_target_seqs for blastn. I'd like to break the query sequences into small segments, and then convert the alignment positions to the origin coordinate system.
#   Many my own scripts will be used.
# [4/21/2022] Use megablast instead of dc-megablast to speed up.
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "bn_task:s", # megablast
);

my $htxt = <<HH;
################################################################################
perl $0 out_prefix 10000 to_rm_contamination.fa db_name


HH

$opts{'bn_task'} //= 'megablast';


!@ARGV and &LogInforSunhh::usage($htxt);

my $opref = shift;
my $segL  = shift;
my $inFa  = shift;
my $dbfn  = shift;

my @torm;

my $stepL = int($segL/2); $stepL < 1 and $stepL = 1;

&runCmd("deal_fasta.pl -chop_seq -chop_len $segL -chop_step $stepL -chop_min 1 $inFa > $opref.seg.fa"); push(@torm, "$opref.seg.fa");
&runCmd("deal_fasta.pl -chop_seq -chop_len $segL -chop_step $stepL -chop_min 1 $inFa -chop_agp > $opref.seg2ori.agp"); push(@torm, "$opref.seg2ori.agp");
&runCmd("deal_fasta.pl -attr key:len $inFa > $opref.ori.kl"); push(@torm, "$opref.ori.kl");

{
  # my $cmd = "blastn -query $opref.seg.fa -out $opref.seg.toDb.bn6 -db $dbfn -task dc-megablast ";
  my $cmd = "blastn -query $opref.seg.fa -out $opref.seg.toDb.bn6 -db $dbfn -task $opts{'bn_task'} ";
  $cmd .= " -evalue 1e-5 -num_threads 80 -max_hsps 20 -max_target_seqs 20  ";
  $cmd .= " -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand staxids sscinames sskingdoms stitle' ";
  &runCmd($cmd); push(@torm, "$opref.seg.toDb.bn6");
  $cmd = "perl /home/Sunhh/tools/github/NGS_data_processing/assemble_tools/cnvt_loc_fromAGP_toAGP_forLoci.pl ";
  $cmd .= " -new_agp $opref.seg2ori.agp  ";
  $cmd .= " -old_loc $opref.seg.toDb.bn6 ";
  $cmd .= " -new_loc $opref.ori.toDb.bn6 ";
  $cmd .= " -colN_seqID 0 ";
  $cmd .= " -colN_seqP 6,7 ";
  &runCmd($cmd); push(@torm, "$opref.ori.toDb.bn6");
  &runCmd("ColLink.pl $opref.ori.toDb.bn6 -f1 $opref.ori.kl -keyC1 0 -keyC2 0 -add -fill NNNN -Col1 1 | deal_table.pl -column 0-11,19,13-18 > $opref.bn6");
}

for (@torm) {
  &fileSunhh::_rmtree($_);
}
#&runCmd("deal_fasta.pl -keep_len 0-$max_len $ctgFa > $opref.tochk.fa");

