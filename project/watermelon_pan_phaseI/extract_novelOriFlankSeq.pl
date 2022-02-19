#!/usr/bin/perl
# [2/7/2022] Fix a bug when the extended sequene of two novel contigs are the same, more than one exactly same sequences will be output.
# [2/8/2022] Although the extended sequence shared by two novel contigs should not be ouput twice in fasta file, its relationship to each novel contig is unique and should be output to ext2nov.tbl.
#   I should add accession ID in the assembled contig ID to distinguish it from those from other accessions.
use strict; 
use warnings;
use LogInforSunhh; 
use fileSunhh;
use fastaSunhh; 
my $fs_obj = fastaSunhh->new();

my $help_txt = <<HH;
perl $0  100000(flank_length)   map.dedup_rmcont_asm.tbl   map.ID1_spec_procID_asmFn_asmPath  out_prefix
# Be aware that the files listed in list "map.ID1_spec_procID_asmFn_asmPath" should exist too!
# Outputs: 
#   out_prefix.ext.fa
#   out_prefix.ori2ext.tbl : In format of deal_fasta.pl -drawLcol ;
#   out_prefix.ext2nov.tbl : In format of deal_fasta.pl -drawLcol ;
HH
!@ARGV and die "$help_txt";

my $ext_len    = shift;
my $f1_nov2asm = shift;
my $f2_asm2fas = shift;
my $opref      = shift;

my $wdir = &fileSunhh::new_tmp_dir( 'create' => 1 );
&fileSunhh::write2file("$wdir/ext.fa", '', '>');
&fileSunhh::write2file("$wdir/ori2ext.tbl", '', '>');
&fileSunhh::write2file("$wdir/ext2nov.tbl", '', '>');

open F1,'<',"$f1_nov2asm" or die ; 
my (%novINasm, %needAccID);
while (<F1>) {
  chomp; 
  my @ta=split(/\t/, $_);
  $ta[0] eq 'dedup_ID' and next;
  $novINasm{$ta[0]} = [ $ta[4], $ta[5], $ta[6], $ta[7] ]; # [acc_ID, asm_ID, asm_start, asm_end]
  $needAccID{$ta[4]} = 1;
}
close F1;
open F2,'<',"$f2_asm2fas" or die;
my %asmSeq;
my %asmID_shrt2long;
while (<F2>) {
  chomp;
  my @ta=split(/\t/, $_);
  defined $ta[2] or next;
  defined $needAccID{$ta[2]} or next;
  $asmSeq{$ta[2]}= $fs_obj->save_seq_to_hash( 'faFile' => $ta[4] );
  for my $k1 (keys %{$asmSeq{$ta[2]}}) {
    $asmSeq{$ta[2]}{$k1}{'seq'} =~ s!\s!!g;
    $asmSeq{$ta[2]}{$k1}{'len'} = length($asmSeq{$ta[2]}{$k1}{'seq'}); 
    if      ( $k1 =~ m!^(NODE_\d+)_length_\d+_cov_\d+\.\d+$! ) {
      # spades_raw: NODE_1303_length_44967_cov_1.754839
      $asmID_shrt2long{$ta[2]}{$1} = $k1;
    } elsif ( $k1 =~ m!^(NODE_\d+)_length_\d+_cov_\d+\.\d+__(\d+)\_(\d+)$! ) {
      # Not very possible.
      $asmID_shrt2long{$ta[2]}{"$1__$2_$3"} = $k1;
    } else {
      $asmID_shrt2long{$ta[2]}{$k1} = $k1;
    }
  }
}
close F2;

# step 3. Extract sequences with flanking region.
my %used_ext;
for my $dedupID (sort keys %novINasm) {
  my ($accID, $asmID, $asmS, $asmE) = @{$novINasm{$dedupID}};
  my $asmID_long = $asmID_shrt2long{$accID}{$asmID};
  defined $asmSeq{$accID}{$asmID_long} or &stopErr("[Err] Failed to find seq for [$dedupID] in [$accID $asmID $asmID_long]\n");
  my ($flank_s, $flank_e) = ($asmS-$ext_len, $asmE+$ext_len);
  $flank_s < 1 and $flank_s = 1;
  $flank_e > $asmSeq{$accID}{$asmID_long}{'len'} and $flank_e = $asmSeq{$accID}{$asmID_long}{'len'};
  my $extID = "${accID}__${asmID}__${flank_s}_${flank_e}_ext";
  &fileSunhh::write2file("$wdir/ext2nov.tbl", join("\t", $extID, $asmS-$flank_s+1, $asmE-$flank_s+1, "+", $dedupID, $accID)."\n", '>>');
  defined $used_ext{$extID} and next;
  $used_ext{$extID} = 1;
  &fileSunhh::write2file("$wdir/ori2ext.tbl", join("\t", "${accID}__$asmID", $flank_s,         $flank_e,         "+", $extID,   $accID)."\n", '>>');
  my $extSeq = substr($asmSeq{$accID}{$asmID_long}{'seq'}, $flank_s-1, $flank_e-$flank_s+1);
  &fileSunhh::write2file("$wdir/ext.fa", ">$extID\n$extSeq\n", '>>');
}

&fileSunhh::_move("$wdir/ext.fa", "${opref}.ext.fa");
&fileSunhh::_move("$wdir/ori2ext.tbl", "${opref}.ori2ext.tbl");
&fileSunhh::_move("$wdir/ext2nov.tbl", "${opref}.ext2nov.tbl");

&fileSunhh::_rmtree($wdir);

# perl /data/Sunhh/watermelon/08.pan_phase_1/map_accDedup_to_asm.pl Clean_UA.CD90.NR.fa.map2Src.tbl > map.dedup_rmcont_asm.tbl
# 
# dedup_ID        rmcont_ID       rmcont_start    rmcont_end      acc_ID  asm_ID  asm_start       asm_end
# Grif_15897_NODE_1303__1-536H    Grif_15897_NODE_1303__1-536H    1       536     Grif_15897      NODE_1303       1       536
# Grif_15897_NODE_18737__1703-2670T       Grif_15897_NODE_18737__1703-2670T       1       968     Grif_15897      NODE_18737      1703    2670
# 
# [Sunhh@panda 08.pan_phase_1]$ cat /data/Sunhh/watermelon/08.pan_phase_1/map.ID1_spec_procID_asmFn_asmPath | awk '$2 == "CA"' | head -4
# ASM1005_PI482246        CA      ASM1005_PI482246        ASM1005_PI482246.ref.fa /data/Sunhh/watermelon/08.pan_phase_1/common_data/ASM1005_PI482246.ref.fa
# WM1060_PI482246 CA      WM1060_PI482246 NA      NA
# WM1191_PI296341 CA      WM1191_PI296341 WM1191_PI296341.spades.500.fa.gz        /data/Sunhh/watermelon/08.pan_phase_1/CA/03_spades_asm/WM1191_PI296341.spades.500.fa.gz
# WM43_PI505586   CA      WM43_PI505586   WM43_PI505586.spades.500.fa.gz  /data/Sunhh/watermelon/08.pan_phase_1/CA/03_spades_asm/WM43_PI505586.spades.500.fa.gz


