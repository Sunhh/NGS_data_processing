#!/usr/bin/perl
# Bedtools and some scripts  are prerequisite.
# [3/24/2022] Count summary information for gene models transferred with Liftoff.
#   R - the genome to transfer from; This is called query in Liftoff.
#   Q - the genome to transfer to; This is called reference in Liftoff.
#   The required counts are:
#     (1) The information of Q gene in Q genome, including
#       old_mrnaID, CDS length, CDS-containing exon number, loc in Q genome;
#     (2) The information of R gene in R genome, including
#       old_mrnaID, CDS length, CDS-containing exon number, loc in R genome;
#     (3) The information of R gene in Q genome, including
#       new_mrnaID, old_mrnaID, CDS length, CDS-containing exon number, liftoff coverage, liftoff identity, loc in Q genome;
#     (4) The mapping relationship from R genes to Q genes (R-to-Q), including
#       R_ID vs Q_ID (. if none), overlapping size (. if no Q gene), 
#   Steps:
#     (1) Get information of Q gene in Q genome and R gene in R genome.
#     (2) Get information of R gene in Q genome after transferring.
#     (3) Get mapping relationship.
#       (3.1) Break R_to_Q.liftoff.gff3 into p.CDS.bed and m.CDS.bed;
#       (3.2) Break Qown.gff3 into p.CDS.bed and m.CDS.bed;
#       (3.3) Run bedtools to get intersection (overlapping) information from R_to_Q.
#       (3.4) Count the features required.
# Although there are steps that use the same input file and can be counted together, I prefer to do it twice to make the functions separate in coding blocks.
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;

!@ARGV and die "perl $0  R_own.gff3  Q_own.gff3  R_to_Q.liftoff.gff3  >  R_to_Q.info.tbl\n";

my ($rOwnGff, $qOwnGff, $r2qGff) = @ARGV;

my $pl_dealGff    = 'perl /home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl';
my $pl_gff2cdsBed = 'perl /data/Sunhh/try_liftoff/tools/cnvt_gff_to_cdsBed.pl';
my $pl_infLOgff   = 'perl /data/Sunhh/try_liftoff/tools/info_liftoffGff.pl';
my $pl_infBTins   = 'perl /data/Sunhh/try_liftoff/tools/info_bedtools_intersect.pl';
my $exe_bedtools  = 'bedtools';

my $wdir = &fileSunhh::new_tmp_dir('create' => 1);
&runCmd("$pl_dealGff -inGff $qOwnGff -getJnLoc > $wdir/qown.gff.JnLoc");
&runCmd("$pl_dealGff -inGff $rOwnGff -getJnLoc > $wdir/rown.gff.JnLoc");

# Step (1) (2):
my (%inf_qOwnG, %inf_rOwnG);
&load_ownG(\%inf_qOwnG, "$wdir/qown.gff.JnLoc");
&load_ownG(\%inf_rOwnG, "$wdir/rown.gff.JnLoc");

# Step (3):
&runCmd("$pl_infLOgff  $r2qGff > $wdir/r2q.liftoff.gff3.tbl");
#     (3) The information of R gene in Q genome, including
#       new_mrnaID, old_mrnaID, CDS length, CDS-containing exon number, liftoff coverage, liftoff identity, loc in Q genome;
my %inf_rInQ;
&load_rInQ(\%inf_rInQ, "$wdir/r2q.liftoff.gff3.tbl");
&runCmd("$pl_gff2cdsBed  $wdir/r2q  $r2qGff ");
&runCmd("$pl_gff2cdsBed  $wdir/qOwn $qOwnGff");
&runCmd("$exe_bedtools intersect -a $wdir/r2q.p.CDS.bed -b $wdir/qOwn.p.CDS.bed -wao > $wdir/is.p");
&runCmd("$exe_bedtools intersect -a $wdir/r2q.m.CDS.bed -b $wdir/qOwn.m.CDS.bed -wao > $wdir/is.m");
# Load bedtools.intersect
&runCmd("$pl_infBTins  $wdir/is.p $wdir/is.m > $wdir/is.b.ovl");
my %ovl_r2q;
&load_r2q(\%ovl_r2q, "$wdir/is.b.ovl");
# Step (3.4): Output summary information.
my @header;
push(@header, qw/trans_R_ID     coverage     identity    R2Q_ovlLen  trans_R_size trans_R_exNum/);
push(@header, qw/ori_R_ID       ori_R_size   ori_R_exNum/);
push(@header, qw/Q_ID           Q_size       Q_exNum/);
push(@header, qw/trans_R_chrIID trans_R_chrS trans_R_chrE trans_R_chrStr/);
push(@header, qw/Q_chrID        Q_chrS       Q_chrE       Q_chrStr/);
push(@header, qw/ori_R_chrID    ori_R_chrS   ori_R_chrE   ori_R_chrStr/);
push(@header, qw/trans_R_CDSblks  Q_CDSblks  ori_R_CDSblks/);
print STDOUT join("\t", @header)."\n";
for my $new_rID (sort { $inf_rInQ{$a}[0] <=> $inf_rInQ{$b}[0] } keys %inf_rInQ) {
  my ($rank_rInQ, $cdsLen_rInQ, $exN_rInQ, $chrID_rInQ, $chrStr_rInQ, $chrS_rInQ, $chrE_rInQ, $cdsblks_rInQ, $coverage, $identity, $ori_rID) = @{$inf_rInQ{$new_rID}};
  # defined $inf_rOwnG{$ori_rID} or die "$ori_rID aa\n";
  my ($rank_rOG,  $cdsLen_rOG,  $exN_rOG,  $chrID_rOG,  $chrStr_rOG,  $chrS_rOG,  $chrE_rOG,  $cdsblks_rOG) = @{$inf_rOwnG{$ori_rID}};
  for my $a1 (@{$ovl_r2q{$new_rID}}) {
    my ($tgt_mrnaID, $trans_cdsLen, $trans_ovlLen) = @$a1;
    my ($rank_qOG, $cdsLen_qOG, $exN_qOG, $chrID_qOG, $chrStr_qOG, $chrS_qOG, $chrE_qOG, $cdsblks_qOG) = (".") x 8;
    if ($tgt_mrnaID ne '.') {
      ($rank_qOG, $cdsLen_qOG, $exN_qOG, $chrID_qOG, $chrStr_qOG, $chrS_qOG, $chrE_qOG, $cdsblks_qOG) = @{$inf_qOwnG{$tgt_mrnaID}};
    }
    my @o1;
    push(@o1, $new_rID,    $coverage,   $identity,  $trans_ovlLen,   $cdsLen_rInQ, $exN_rInQ);
    push(@o1, $ori_rID,    $cdsLen_rOG, $exN_rOG);
    push(@o1, $tgt_mrnaID, $cdsLen_qOG, $exN_qOG);
    push(@o1, $chrID_rInQ, $chrS_rInQ,  $chrE_rInQ, $chrStr_rInQ);
    push(@o1, $chrID_qOG,  $chrS_qOG,   $chrE_qOG,  $chrStr_qOG);
    push(@o1, $chrID_rOG,  $chrS_rOG,   $chrE_rOG,  $chrStr_rOG);
    push(@o1, $cdsblks_rInQ, $cdsblks_qOG, $cdsblks_rOG);
    print STDOUT join("\t", @o1)."\n";
  }
}
&fileSunhh::_rmtree($wdir);

sub load_r2q {
  my ($hR, $inF) = @_;
  my $cnt_tmp = 0;
  for my $l1 (&fileSunhh::load_tabFile($inF)) {
    $cnt_tmp ++;
    $l1->[0] eq 'trans_mrnaID' and next;
    push(@{$hR->{$l1->[0]}}, [@{$l1}[1,2,3]]); # [tgt_mrnaID, trans_cdsLen, trans_ovlLen]
  }
}# load_r2q()

sub load_rInQ {
  my ($hR, $inF) = @_;
  my $cnt_tmp = 0;
  for my $l1 (&fileSunhh::load_tabFile($inF)) {
    $cnt_tmp ++;
    $l1->[0] eq 'new_mrnaID' and next;
    $hR->{$l1->[0]} = [$cnt_tmp, @{$l1}[3, 4,  8,11,9,10,  12, 1, 2, 5]  ]; # [rank_cnt, cds_len, CDSBlocksNum, chrID,strand,chrS,chrE, CDSBlocks, coverage, identity, old_mrnaID]
  }
}# load_rInQ()
sub load_ownG {
  my ($hR, $inF) = @_;
  my $cnt_tmp = 0;
  for my $l1 (&fileSunhh::load_tabFile($inF)) {
    $cnt_tmp ++;
    $l1->[0] eq 'mrnaID' and next;
    $hR->{$l1->[0]} = [$cnt_tmp, $l1->[8], $l1->[10], @{$l1}[2,5,6,7], $l1->[9]]; # [rank_cnt, cds_len, CDSBlocksNum, chrID,strand,chrS,chrE, CDSBlocks]
  }
}# load_ownG()


