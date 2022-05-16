#!/usr/bin/perl -w
# [5/16/2022] Fix a bug which causes the missing of some genes that are not in a selected syntenic block but grouped with genes in a syntenic block.
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;
use mathSunhh;
my $ms_obj = mathSunhh->new();
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "in_grp_mat:s",   # ov2/comb.grp2.novl_loc: Grp_name_0 / size_1 / tag1:gene1 / tag1:gene2 / tag2:gene1 / tag3:gene1 / ...
    # The tag should be two-letter long.
  "in_geneJnLoc:s", # input/all.gff3JnLoc:    mrnaID_0   / geneID / SeqID_2 / mrnaStart / mrnaEnd / mrnaStrand_5 / CDSStart_6 / CDSEnd_7 / LenInCDS / ...
  "in_orthLis:s",   # input/all.syn.tbl.orth: BlkID      / Chrom1 / Start1  / End1      / Chrom2  / Start2       / End2       / Strand   / AlnScore / AlnEvalue / AlnNumber
  "in_novGenLis:s",
  "opref:s",
  "min_anyOvlR:f",
  "help!"
);


$opts{'opref'} //= 'outGrp';
my $default_min_anyOvlR = 0.5;

my $htxt = <<HH;
################################################################################
perl $0  -in_grp_mat ov2/comb.grp2.novl_loc  -in_geneJnLoc input/all.gff3JnLoc  -in_orthLis input/all.syn.tbl.orth -opref out

 -in_grp_mat     [file path] the output file from remove_ovl_loc.pl.
 -in_geneJnLoc   [file path] combined .gff3.JnLoc
 -in_orthLis     [file path] combined mcscan.collinearity.tbl.orth
 -opref          [output prefix]

 -min_anyOvlR    [0-1] Default $default_min_anyOvlR . Minimum ratio of overlapping region to the shorter block length.

 -in_novGenLis   [file path] A file with the first column as the mrnaIDs of novel genes in the pan-genome.
                   Currently, I'll keep all the novel genes in the first synteny block, which may not be so accurate if a novel gene is brought in by other genes.

HH

$opts{'min_anyOvlR'} //= $default_min_anyOvlR;
$opts{'min_anyOvlR'} > 0 or &stopErr("[Err] -min_anyOvlR must be bigger than zero.\n");

for my $required (qw/in_grp_mat in_geneJnLoc in_orthLis/) {
  defined $opts{$required} or &LogInforSunhh::usage($htxt);
}

# Load information.
my %ha_gLoc = &load_geneJnLoc($opts{'in_geneJnLoc'}); # {mrnaID} = [SeqID, CDSStart, CDSEnd]
my %ha_sBlk = &load_orthLis($opts{'in_orthLis'});     # {chrID_1}{chrID_2} = [newBlkID, Start1, End1, Start2, End2]
my @ar_grp = &load_grpMat($opts{'in_grp_mat'});       # ([ [tag1, gene1], [tag2, gene2], [tag2, gene3], ...], ...)
my %ha_novGen = ();
defined $opts{'in_novGenLis'} and %ha_novGen = &load_novGenLis($opts{'in_novGenLis'}); # {mrnaID} = 1;

# Revise groups one by one.
my @connected_grp; # Groups of genes that can be connected by syntenic blocks.
my @sep_grp; # Groups of genes that are not in any syntenic blocks.
my @orphan_grp;
for my $gt (@ar_grp) {
  scalar(@$gt) == 1 and do { push(@orphan_grp, [@$gt]); next; };
  my @old_gt = @$gt;
  my @new_gt;
  my %blk_pairs;
  my @blk_pkey;
  my @in_blks;
  my %in_genesH;
  my %in_pkey;
  my %rest_genesH;
  my @rest_gt;
  REVISE_GRP:
  # Count the number of gene pairs carried in each synteny block.
  #   In normal case, each synteny block should be supported only once,
  #   Usually the block with high gene pairs should contain tandem duplicate genes, and I want to select this block as the seed.
  for (my $i=0; $i<@old_gt; $i++) {
    my $tag1  = $old_gt[$i][0];
    my $gene1 = $old_gt[$i][1];
    if (defined $ha_novGen{$gene1}) {
      $in_genesH{"$tag1\t$gene1"} = 1;
      next;
    }
    my ($chr1, $s1, $e1) = &ret_gene_info($gene1, \%ha_gLoc); # ( chrID, start, end )
    defined $ha_sBlk{$chr1} or next;
    for (my $j=$i+1; $j<@old_gt; $j++) {
      my $tag2  = $old_gt[$j][0];
      my $gene2 = $old_gt[$j][1];
      my ($chr2, $s2, $e2) = &ret_gene_info($gene2, \%ha_gLoc);
      defined $ha_sBlk{$chr1}{$chr2} or next;
      for my $tt1 (@{$ha_sBlk{$chr1}{$chr2}}) {
        ( $tt1->[1] < $s1 and $tt1->[2] > $e1 ) or next; # Fully included, ignoring those mapped locations spanning the boundaries.
        ( $tt1->[3] < $s2 and $tt1->[4] > $e2 ) or next;
        my $pkey = ($chr1 cmp $chr2) ? "$chr1\t$chr2\t$tt1->[0]" : "$chr2\t$chr1\t$tt1->[0]";
        push(@{$blk_pairs{$pkey}}, [$gene1, $gene2, $tag1, $tag2]);
      }
    }
  }
  # Build a new group with blocks that can be connected to the seed block.
  @blk_pkey = sort { scalar(@{$blk_pairs{$b}}) <=> scalar(@{$blk_pairs{$a}}) || $a cmp $b } keys %blk_pairs;
  if (scalar(@blk_pkey) == 0) {
    push(@sep_grp, [@old_gt]);
    next;
  }
  @in_blks = ( [&cnvt_pkey2arr($blk_pkey[0])] ); # ( [chr1, chr2, blkID], ...)
  my $prevN_in_blks = -1;
  my $afterN_in_blks = -2;
  my %ha_in_blks = ( 0 => 1 );
  while ($prevN_in_blks != $afterN_in_blks) {
    $prevN_in_blks = scalar(@in_blks);
    for (my $i=1; $i<@blk_pkey; $i++) {
      defined $ha_in_blks{$i} and next;
      my ($chk_chr1, $chk_chr2, $chk_blkID) = &cnvt_pkey2arr($blk_pkey[$i]);
      my $is_ovlBlk = &if_ovlBlk( \@in_blks, [ $chk_chr1, $chk_chr2, $chk_blkID ], \%ha_sBlk, $opts{'min_anyOvlR'} );
      $is_ovlBlk == 1 and do { push(@in_blks, [ $chk_chr1, $chk_chr2, $chk_blkID ] ); $ha_in_blks{$i} = 1; };
    }
    $afterN_in_blks = scalar(@in_blks);
  }
  # Put genes in the new group blocks into a new gene group.
  for my $bt (@in_blks) {
    my $pkey = "$bt->[0]\t$bt->[1]\t$bt->[2]";
    $in_pkey{$pkey} = 1;
    defined $blk_pairs{$pkey} or &stopErr("[Err] No data for pkey [$pkey]\n");
    for my $pt (@{$blk_pairs{$pkey}}) {
      $in_genesH{"$pt->[2]\t$pt->[0]"} = 1;
      $in_genesH{"$pt->[3]\t$pt->[1]"} = 1;
    }
  }
  @new_gt = map { [split("\t", $_)] } sort keys %in_genesH;
  # Put the rest genes into another new group for another check, and add the novel genes to in_genesH group.
  for my $a1 (@old_gt) {
    my $a2 = "$a1->[0]\t$a1->[1]";
    defined $in_genesH{$a2} or $rest_genesH{$a2} = 1;
  }
  @rest_gt = map { [split("\t", $_)] } sort keys %rest_genesH;
  # Do this division again if there is a new gene group found.
  if (scalar(@new_gt) > 0) {
    push(@connected_grp, [@new_gt]);
    if (scalar(@rest_gt) > 0) {
      # Initialize.
      @old_gt = @rest_gt;
      @new_gt = %blk_pairs = @blk_pkey = @in_blks = %in_genesH = %in_pkey = %rest_genesH = @rest_gt = ();
      goto REVISE_GRP;
    }
  } else {
    scalar(@rest_gt) > 0 or &stopErr("[Err] Impossible 1:\n");
    push(@sep_grp, [@rest_gt]);
  }
}

# Output results.
### ([ [tag1, gene1], [tag2, gene2], [tag2, gene3], ...], ...)
{
  my $gCnt = 0;
  @connected_grp = sort { scalar(@$b) <=> scalar(@$a) } @connected_grp;
  @sep_grp = sort { scalar(@$b) <=> scalar(@$a) } @sep_grp;
  my $gN = scalar(@connected_grp) + scalar(@sep_grp);
  my $gCntLen = length($gN)+1;
  my $ofh1 = &openFH("$opts{'opref'}.grp.tbl", '>');
  for my $pt (@connected_grp) {
    $gCnt ++;
    my $grpID = sprintf("GrpSyn_%0${gCntLen}d", $gCnt);
    print ${ofh1} join("\t", $grpID, scalar(@$pt), map { "$_->[0]:$_->[1]" } @$pt)."\n";
  }
  for my $pt (@sep_grp) {
    $gCnt ++;
    my $grpID = sprintf("GrpSep_%0${gCntLen}d", $gCnt);
    print ${ofh1} join("\t", $grpID, scalar(@$pt), map { "$_->[0]:$_->[1]" } @$pt)."\n";
  }
  for my $pt (@orphan_grp) {
    $gCnt ++;
    my $grpID = sprintf("GrpOrp_%0${gCntLen}d", $gCnt);
    print ${ofh1} join("\t", $grpID, scalar(@$pt), map { "$_->[0]:$_->[1]" } @$pt)."\n";
  }
}
# All done.

####################################################################################################
# Sub-routines.
sub if_ovlBlk {
  my ($poolBlk, $testBlk, $blkHR, $min_anyOvlR) = @_;
  my $is_ovl = 0;
  my ($ts1, $te1, $ts2, $te2) = &ret_se_byChrBlkID($testBlk->[0], $testBlk->[1], $testBlk->[2], $blkHR);
  for my $tp (@$poolBlk) {
    my ($ps1, $pe1, $ps2, $pe2) = &ret_se_byChrBlkID($tp->[0], $tp->[1], $tp->[2], $blkHR);
    if ($testBlk->[0] eq $tp->[0]) {
      my ($ovlLen) = $ms_obj->ovl_region($ts1, $te1, $ps1, $pe1);
      ($ovlLen >= $min_anyOvlR * ($te1-$ts1+1) or $ovlLen >= $min_anyOvlR * ($pe1-$ps1+1) ) and do {$is_ovl = 1; last;};
    }
    if ($testBlk->[0] eq $tp->[1]) {
      my ($ovlLen) = $ms_obj->ovl_region($ts1, $te1, $ps2, $pe2);
      ($ovlLen >= $min_anyOvlR * ($te1-$ts1+1) or $ovlLen >= $min_anyOvlR * ($pe2-$ps2+1) ) and do {$is_ovl = 1; last;};
    }
    if ($testBlk->[1] eq $tp->[0]) {
      my ($ovlLen) = $ms_obj->ovl_region($ts2, $te2, $ps1, $pe1);
      ($ovlLen >= $min_anyOvlR * ($te2-$ts2+1) or $ovlLen >= $min_anyOvlR * ($pe1-$ps1+1) ) and do {$is_ovl = 1; last;};
    }
    if ($testBlk->[1] eq $tp->[1]) {
      my ($ovlLen) = $ms_obj->ovl_region($ts2, $te2, $ps2, $pe2);
      ($ovlLen >= $min_anyOvlR * ($te2-$ts2+1) or $ovlLen >= $min_anyOvlR * ($pe2-$ps2+1) ) and do {$is_ovl = 1; last;};
    }
  }
  return($is_ovl);
}# if_ovlBlk()

sub ret_se_byChrBlkID {
  my ($chr1, $chr2, $blkID, $blkHR) = @_;
  for my $t1 (@{$blkHR->{$chr1}{$chr2}}) {
    $t1->[0] == $blkID or next;
    return(@{$t1}[1,2,3,4]);
  }
}# ret_se_byChrBlkID()

sub cnvt_pkey2arr {
  my ($str) = @_;
  $str =~ m!^(\S+)\t(\S+)\t(\d+)$! or &stopErr("[Err] Bad pkey [$str]\n");
  return($1,$2,$3);
}# cnvt_pkey2arr()

# Return: ( chrID, start, end )
sub ret_gene_info {
  my ($geneID, $glocH) = @_;
  if (defined $glocH->{$geneID}) {
    return( $glocH->{$geneID}[0], $glocH->{$geneID}[1], $glocH->{$geneID}[2] );
  } elsif ($geneID =~ m!^(\S+):(\d+)\-(\d+):[+-]$!) {
    # CicolChr04:31896151-31897083:-
    return( $1, $2, $3 );
  } else {
    &stopErr("[Err] Failed to parse gene ID [$geneID]\n");
  }
}# ret_gene_info()

#  "in_orthLis:s",   # input/all.syn.tbl.orth: BlkID      / Chrom1 / Start1  / End1      / Chrom2  / Start2       / End2       / Strand   / AlnScore / AlnEvalue / AlnNumber
# Return: %back;
#   {chrID_1}{chrID_2} = [newBlkID, Start1, End1, Start2, End2]
#   {chrID_2}{chrID_1} = [newBlkID, Start2, End2, Start1, End1]
sub load_orthLis {
  my ($fn) = @_;
  my %back;
  my $fh = &openFH($fn, '<');
  while (<$fh>) {
    chomp;
    my @ta=split(/\t/, $_);
    $ta[0] eq 'BlkID' and next;
    push(@{$back{$ta[1]}{$ta[4]}}, [$., @ta[2,3,5,6]]);
    if ($ta[1] ne $ta[4]) {
      push(@{$back{$ta[4]}{$ta[1]}}, [$., @ta[5,6,2,3]]);
    }
  }
  close($fh);
  return(%back);
}# load_orthLis()

#  "in_geneJnLoc:s", # input/all.gff3JnLoc:    mrnaID_0   / geneID / SeqID_2 / mrnaStart / mrnaEnd / mrnaStrand_5 / CDSStart_6 / CDSEnd_7 / LenInCDS / ...
# Return: %bb;
#   {mrnaID} = [SeqID, CDSStart, CDSEnd]
sub load_geneJnLoc {
  my ($fn) = @_;
  my %back;
  my $fh = &openFH($fn, '<');
  while (<$fh>) {
    chomp;
    my @ta=split(/\t/, $_);
    $ta[0] eq 'mrnaID' and next;
    defined $back{$ta[0]} and &stopErr("[Err] Repeat mrnaID [$ta[0]]\n");
    $back{$ta[0]} = [@ta[2,6,7]];
  }
  close($fh);
  return(%back);
}# load_geneJnLoc()

#  "in_grp_mat:s",   # ov2/comb.grp2.novl_loc: Grp_name_0 / size_1 / tag1:gene1 / tag1:gene2 / tag2:gene1 / tag3:gene1 / ...
# Return: @bb;
#   ([tag1, gene1], [tag2, gene2], [tag2, gene3], ...)
sub load_grpMat {
  my ($fn) = @_;
  my @back; # {grpID} = [ [tag1, gene1], [tag2, gene2], ... ]
  my $fh = &openFH($fn, '<');
  while (<$fh>) {
    chomp;
    my @ta=split(/\t/, $_);
    $ta[0] eq 'Grp_name' and next;
    my @tb;
    for my $b1 (@ta[2 .. $#ta]) {
      $b1 =~ m!^([^\s:]+):(\S+)$! or &stopErr("[Err] bad format in grp_mat: |$b1|\n");
      push(@tb, [$1, $2]);
    }
    scalar(@tb) > 0 or next;
    push(@back, [@tb]);
  }
  close($fh);
  return(@back);
}# load_grpMat()

sub load_novGenLis {
  my ($fn) = @_;
  my $fh = &openFH($fn, '<');
  my %back;
  while (<$fh>) {
    chomp;
    my @ta=&splitL("\t", $_);
    $back{$ta[0]} = 1;
  }
  close($fh);
  return(%back);
}# load_novGenLis()

