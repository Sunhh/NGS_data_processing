#!/usr/bin/perl
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "gene_distance:i", # 2000
  "help!"
);

$opts{'gene_distance'} //= 2000;

my $htxt = <<HH;
################################################################################
# perl $0   out_prefix   input.gff3   input.p.fas   db_diamond
# 
# ## Use this script to retrieve gene models that
#   (1) have complete protein sequences in db_diamond;
#   (2) do not overlap each other with at least overall 80% identity or overall 50% sum-up matches 
#       (details determined by perl script rmRedunt_inputProt.pl).
#   (3) separate at least 'gene_distance' ($opts{'gene_distance'} bp) away from each other.
#
################################################################################
HH

scalar(@ARGV) == 4 or &LogInforSunhh::usage($htxt);

my $opref = shift;
my $fn_inGff = shift;
my $fn_inFas = shift;
my $db       = shift;

# Required paths
my $exe_diamond = 'diamond';
my $pl_dealFas  = 'deal_fasta.pl';
my $pl_dealGff  = 'perl /home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl';
my $pl_rmRedPr  = 'perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/rmRedunt_inputProt.pl';


my $wdir = &fileSunhh::new_tmp_dir('create' => 1);

{
  my $cmd = "$exe_diamond blastp ";
  $cmd .= " -q $fn_inFas ";
  $cmd .= " -d $db ";
  $cmd .= " --masking none ";
  $cmd .= " -e 1e-4 -k 20 -p 20 ";
  $cmd .= " -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen ";
  $cmd .= " -o $wdir/bp2db.bp6 ";
  &runCmd($cmd);
  $cmd = '';
  open F1,'<',"$wdir/bp2db.bp6" or die "__1__\n";
  open O1,'>',"$wdir/bp2db.bp6.comID" or die "__1b__\n";
  my %h;
  while (<F1>) {
    chomp;
    my @ta=split(/\t/, $_);
    defined $h{$ta[0]} and next;
    $ta[6] <= 5 or next;
    $ta[7] >= $ta[12]-5 or next;
    $ta[7]-$ta[6]+1 >= 0.95 * $ta[12] or next;
    $ta[8] <= 5 or next;
    $ta[9] >= $ta[13]-5 or next;
    $ta[9]-$ta[8]+1 >= 0.95 * $ta[13] or next;
    $h{$ta[0]} = 1;
    print O1 "$ta[0]\n";
  }
  close O1;
  close F1;
  &runCmd("$pl_dealFas $fn_inFas -drawByList -drawWhole -drawLcol 0 -drawList $wdir/bp2db.bp6.comID > $wdir/bp2db.bp6.comID.p.fa");
}

{
  my $cmd = "$pl_rmRedPr ";
  $cmd .= " -prot_qry $wdir/bp2db.bp6.comID.p.fa ";
  $cmd .= " -opref    $wdir/rmRed ";
  $cmd .= " -minIdentity 80 -blastp_para ' -evalue 1e-5 -seg no -num_threads 20 -max_target_seqs 30 ' ";
  &runCmd($cmd);
  $cmd = "";
  &runCmd("$pl_dealFas $wdir/bp2db.bp6.comID.p.fa -attr key | tail -n +2 | deal_table.pl -kSrch_drop -kSrch_idxCol 0 -kSrch_srcCol 0 -kSrch_idx $wdir/rmRed.bp6.redund_list > $wdir/rmRed.kept_list");
}

&runCmd("$pl_dealGff -inGff $fn_inGff -gffret $wdir/rmRed.kept_list > $wdir/rmRed.kept.gff3");
&runCmd("$pl_dealGff -inGff $wdir/rmRed.kept.gff3 -islandGene $opts{'gene_distance'} -islandStrand Both > $opref.CompleteRmRedIsland.gff3");
&runCmd("cp -p $wdir/bp2db.bp6.comID $opref.completeProt.ID");
&runCmd("cp -p $wdir/rmRed.kept_list $opref.completeRmRed.ID");
&runCmd("$pl_dealGff -inGff $opref.CompleteRmRedIsland.gff3 -getJnLoc > $opref.CompleteRmRedIsland.gff3.JnLoc");
&runCmd("$pl_dealFas $fn_inFas -drawByList -drawLcol 0 -drawWhole -drawList $opref.CompleteRmRedIsland.gff3.JnLoc > $opref.CompleteRmRedIsland.p.fa");


&fileSunhh::_rmtree($wdir);

