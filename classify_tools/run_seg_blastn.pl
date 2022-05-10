#!/usr/bin/perl
# [4/20/2022] Because the input contigs are super long, it may lower the sensitivity by setting lower -max_hsp and -max_target_seqs for blastn. I'd like to break the query sequences into small segments, and then convert the alignment positions to the origin coordinate system.
#   Many my own scripts will be used.
# [4/21/2022] Use megablast instead of dc-megablast to speed up.
# [5/5/2022]
# [5/10/2022] Fix a bug for score_d2top filtering.
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "bn_task:s", # megablast
  "score_d2top:f", # 0.1
);

my $htxt = <<HH;
################################################################################
perl $0 out_prefix 10000 to_rm_contamination.fa db_name

-score_d2top     [0.1] [0-1)

HH

$opts{'bn_task'}     //= 'megablast';
$opts{'score_d2top'} //= 0.1;

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
  my $cmd = '';
  if ($opts{'bn_task'} eq 'blastx') {
    $cmd = "diamond blastx --masking none -e 1e-5 -k 50 --max-hsps 50 -p 80 ";
    $cmd .= " -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qstrand staxids sscinames sskingdoms stitle ";
    $cmd .= " -q $opref.seg.fa -d $dbfn -o $opref.seg.toDb.bn6"; # I keep the file name although it should be .bx6.
  } else {
    $cmd = "blastn -query $opref.seg.fa -out $opref.seg.toDb.bn6 -db $dbfn -task $opts{'bn_task'} ";
    $cmd .= " -evalue 1e-5 -num_threads 80 -max_hsps 20 -max_target_seqs 20  ";
    $cmd .= " -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand staxids sscinames sskingdoms stitle' ";
  }
  &runCmd($cmd); push(@torm, "$opref.seg.toDb.bn6");
  if ($opts{'bn_task'} eq 'blastx' and $opts{'score_d2top'} < 1 and $opts{'score_d2top'} >= 0) {
    open F1,'<', "$opref.seg.toDb.bn6" or die;
    my (%hits, @qID);
    while (<F1>) {
      chomp;
      my @ta=split(/\t/, $_);
      defined $hits{$ta[0]} or push(@qID, $ta[0]);
      if ($ta[15] eq '' or $ta[16] eq 'N/A') {
        $hits{$ta[0]}{'topScore'} //= $ta[11]; 
        $hits{$ta[0]}{'topScore'} < $ta[11] and $hits{$ta[0]}{'topScore'} = $ta[11];
      } else {
        $hits{$ta[0]}{'topScore'} //= -1;
      }
      $ta[11] >= $hits{$ta[0]}{'topScore'} * (1-$opts{'score_d2top'}) or next;
      push(@{$hits{$ta[0]}{'aln'}}, [@ta]);
    }
    close F1;
    open O1,'>', "$opref.seg.toDb.bn6.1" or die;
    for my $qid (@qID) {
      for (@{$hits{$qid}{'aln'}}) {
        $_->[11] >= $hits{$qid}{'topScore'} * (1-$opts{'score_d2top'}) or next;
        print O1 join("\t", @$_)."\n";
      }
    }
    close O1;
    &fileSunhh::_move("$opref.seg.toDb.bn6.1", "$opref.seg.toDb.bn6");
  }
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

