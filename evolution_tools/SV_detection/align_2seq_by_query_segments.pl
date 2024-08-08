#!/usr/bin/perl
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;

!@ARGV and die "perl $0 query.fa subject.fa out_prefix 'chop_parameters' '-max_diffR 0.01'\n  chop_parameters: -chop_seq -chop_len 1000 -chop_step 500 -chop_min 1000\n";

my $pl_dealFasta='deal_fasta.pl';
my $pl_topBn='perl /home/Sunhh/tools/github/NGS_data_processing/evolution_tools/ortho_tools/filter_bp6_byTopScore.pl';
my $pl_cnvtLoc='perl /home/Sunhh/tools/github/NGS_data_processing/assemble_tools/cnvt_loc_fromAGP_toAGP_forLoci.pl';
my $bn6 = 'blastn -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand"';
$bn6 .= " -num_threads 30 -evalue 1e-5";

my $qFa = shift;
my $sFa = shift;
my $opref = shift;
my $paraChop = '-chop_seq -chop_len 1000 -chop_step 500 -chop_min 1000';
scalar(@ARGV) > 0 and $paraChop = shift;
my $paraTop  = '-max_diffR 0.01';
scalar(@ARGV) > 0 and $paraTop  = shift;

# Prepare database
my $db = $sFa;
my $wrkDir = &fileSunhh::new_tmp_dir('create'=>1);;
unless (-e "$db.nsq") {
  &runCmd("cp -p $sFa $wrkDir/r.fa");
  $db = "$wrkDir/r.fa";
  &runCmd("makeblastdb -in $wrkDir/r.fa -dbtype nucl");
}

# Chop query sequences.
&runCmd("$pl_dealFasta $paraChop $qFa > $wrkDir/qseg.fa");
&runCmd("$pl_dealFasta $paraChop $qFa -chop_agp > $wrkDir/qseg.agp");

# Run blastn
&runCmd("$bn6 -db $db -query $wrkDir/qseg.fa -out $opref.seg.bn6");
&runCmd("$pl_topBn $paraTop $opref.seg.bn6 > $opref.seg.bn6.s1");
&runCmd("$pl_cnvtLoc -new_agp $wrkDir/qseg.agp -old_loc $opref.seg.bn6.s1 -new_loc $opref.s1.bn6 -colN_seqID 0 -colN_seqP 6,7,12 -colN_seqStr 14");

# Done.
&tsmsg("[Red] done. output files are $opref.s1.bn6 and others ($opref.seg.bn6*)\n");

