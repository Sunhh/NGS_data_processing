#!/usr/bin/perl
# 20250114: update to be more flexible.
use strict;
use warnings;
use Getopt::Long;
use SeqAlnSunhh;
use LogInforSunhh;

my %opts;
GetOptions(\%opts,
  "in_bam:s",                # inBam
  "target_loc:s",           # duplicate:66251-66252
  "positive_control_loc:s", # duplicate:131502-131503
  "flank_len:i",             # 2000;
  "max_insertion_len:i",     # 2000;
  "help!",
);
# >duplicate [seq_source=CLV01_Chr01:607964-740465] [break_point=739465-739466]

my $help_txt = <<HH;
######################################################################
perl $0 -in_bam in.bam -target_loc duplicate:66241-66262  -positive_control_loc duplicate:131502-131503

-help
-flank_len         [2000]

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt);
for my $k (qw/target_loc positive_control_loc in_bam/) {
  defined $opts{$k} or &LogInforSunhh::usage($help_txt);
}
$opts{'flank_len'} //= 2000;
$opts{'max_insertion_len'} //= 2000;

my ($tgt_id, $tgt_s, $tgt_e) = ($opts{'target_loc'} =~ m!^\s*(\S+):(\d+)\-(\d+)\s*$!) or die "bad target_loc\n";
my ($pos_id, $pos_s, $pos_e) = ($opts{'positive_control_loc'} =~ m!^\s*(\S+):(\d+)\-(\d+)\s*$!) or die "bad target_loc\n";

my $tgt_seg_s = $tgt_s - $opts{'flank_len'}; $tgt_seg_s < 1 and $tgt_seg_s = 1;
my $tgt_seg_e = $tgt_e + $opts{'flank_len'};
my $pos_seg_s = $pos_s - $opts{'flank_len'}; $pos_seg_s < 1 and $pos_seg_s = 1;
my $pos_seg_e = $pos_e + $opts{'flank_len'};


my %flag_goodR = %{ &SeqAlnSunhh::mk_flag("keep"=>"0=1,2=0,3=0,4=1,5=0") };
my $fn = $opts{'in_bam'};

my $tgt_cnt = 0;
my $pos_cnt = 0;

open F1,"-|", "samtools view $fn $tgt_id:$tgt_seg_s-$tgt_seg_e" or die "$!\n";
while (<F1>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[6] eq "=" or next;
  defined $flag_goodR{$ta[1]} or next;
  my $rP_s = $ta[7];
  my ($rd_len, $ref_span) = &SeqAlnSunhh::cigar_array2len( &SeqAlnSunhh::cigar_str2array( $ta[5] ) );
  my $rP_e = $ta[3] + $ref_span - 1;
  $rP_e-$rP_s+1 <= $opts{'max_insertion_len'} or next;
  if ($rP_s <= $tgt_s and $rP_e >= $tgt_e) {
    $tgt_cnt ++;
  }
}
close F1;

open F2,"-|", "samtools view $fn $pos_id:$pos_seg_s-$pos_seg_e " or die "$!\n";
while (<F2>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[6] eq "=" or next;
  defined $flag_goodR{$ta[1]} or next;
  my $rP_s = $ta[7];
  my ($rd_len, $ref_span) = &SeqAlnSunhh::cigar_array2len( &SeqAlnSunhh::cigar_str2array( $ta[5] ) );
  my $rP_e = $ta[3] + $ref_span - 1;
  $rP_e-$rP_s+1 <= $opts{'max_insertion_len'} or next;
  if ($rP_s <= $pos_s and $rP_e >= $pos_e) {
    $pos_cnt ++;
  }
}
close F2;

print STDOUT join("\t", $fn, $tgt_cnt, $pos_cnt)."\n";

