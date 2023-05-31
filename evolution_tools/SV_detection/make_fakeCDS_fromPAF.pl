#!/usr/bin/perl
# 5/2/2023: Generate 'conserved' anchors for use in anchorwave.
# 5/30/2023: Fix a bug to get unique blocks on query.
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;
use fastaSunhh;
use mathSunhh;

!@ARGV and die "perl $0 out_prefix refU43_22CEXU11.paf refU43.fa\n";

my $anchorLen = 3000; # Anchor length.
my $minAlnLen = 10e3; # Minimum length of a block which can be accepted as a good alignment.
my $maxOvlLen = 1e3;  # Two blocks with overlapping length >= $maxOvlLen is treated as overlapping.

my $opre   = shift;
my $alnPaf = shift;
my $refFa  = shift;

my $ms_obj = mathSunhh->new();
my $fas_obj = fastaSunhh->new();
my %seqR = %{$fas_obj->save_seq_to_hash('faFile'=>$refFa)};
for (keys %seqR) {$seqR{$_}{'seq'} =~ s!\s!!g; $seqR{$_}{'seq'}=uc($seqR{$_}{'seq'}); $seqR{$_}{'len'}=length($seqR{$_}{'seq'});}

# my $wrk_dir = &fileSunhh::new_tmp_dir('create'=>1);

# Load alignments.
my (@aln1, %loc1R, %loc1Q);
&tsmsg("[Msg] Loading PAF alignments.\n");
for my $l1 (&fileSunhh::load_tabFile($alnPaf)) {
  $l1->[11] > 20 or next;
  $l1->[10] > $minAlnLen or next;
  $l1->[2] ++;
  $l1->[7] ++;
  push(@aln1, [@{$l1}[5,7,8, 0,2,3]]);
  push(@{$loc1R{$l1->[5]}}, [@{$l1}[7,8]]);
  push(@{$loc1Q{$l1->[0]}}, [@{$l1}[2,3]]);
}
for (keys %loc1R) {@{$loc1R{$_}} = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @{$loc1R{$_}};}
for (keys %loc1Q) {@{$loc1Q{$_}} = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @{$loc1Q{$_}};}

# Get unique alignments.
&tsmsg(join("", "[Msg] Get unique hits from ", scalar(@aln1), " alignments\n"));
my @aln2;
my %tmp1;
for my $a1 (@aln1) {
  $tmp1{'n_aln1'} ++;
  $tmp1{'n_aln1'} % 100 == 1 and &tsmsg("[Msg]  Processing $tmp1{'n_aln1'} -th alignment.\n");
  my $stime = 0;
  for my $a2 (@{$loc1R{$a1->[0]}}) {
    $a1->[1] > $a2->[1] and next;
    $a1->[2] < $a2->[0] and last;
    my ($ovlLen, $ovlAH) = $ms_obj->ovl_region($a1->[1], $a1->[2], $a2->[0], $a2->[1]);
    $ovlLen >= $maxOvlLen or next;
    $stime ++;
    $stime >= 2 and last;
  }
  $stime >= 2 and next;
  my $qtime = 0;
  for my $a2 (@{$loc1Q{$a1->[3]}}) {
    $a1->[4] > $a2->[1] and next;
    $a1->[5] < $a2->[0] and last;
    my ($ovlLen, $ovlAH) = $ms_obj->ovl_region($a1->[4], $a1->[5], $a2->[0], $a2->[1]);
    $ovlLen >= $maxOvlLen or next;
    $qtime ++;
    $qtime >= 2 and last;
  }
  $qtime >= 2 and next;
  push(@aln2, [@$a1]);
}
@aln2 = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @aln2;
&tsmsg(join('', "[Msg] Got ", scalar(@aln2), " unique hits.\n"));

# Find good anchors in @aln2. And then output files.
### @aln1 is not used anymore.
### To reduce memory usage, stop removing repeat anchors. This step will be done by AnchorWave.
my $prevEnd = -1;
&tsmsg("[Msg] Looking for good anchors.\n");
%tmp1 = ();
my $ofh2=&openFH("$opre.anchor.gff3", '>');
for my $a0 (@aln2) {
  my ($s0, $e0) = ($a0->[1], $a0->[2]);
  $tmp1{'n_aln1'} ++;
  $tmp1{'n_aln1'} % 10 == 1 and &tsmsg("[Msg]  Processing $tmp1{'n_aln1'} -th hit: [@$a0].\n");
  # &tsmsg("[Msg]   Processing $tmp1{'n_aln1'} -th hit: [@$a0].\n");
  $s0 <= $prevEnd and $s0 = $prevEnd+1;
  for (my $s1=$s0; $s1+$anchorLen-1 <= $e0-$anchorLen; $s1+=$anchorLen) {
    my $e1 = $s1+$anchorLen-1;
    $e1+$anchorLen > $e0 and $e1 = $e0;
    my $seq1 = substr($seqR{$a0->[0]}{'seq'}, $s1-1, $e1-$s1+1);
    $seq1 =~ m![nN]! and next;
    my $k=join('_', $a0->[0], $s1, $e1);
    print {$ofh2} join("\t", $a0->[0], ".", "gene", $s1, $e1, qw/. + ./, "ID=${k}_gene")."\n";
    print {$ofh2} join("\t", $a0->[0], ".", "mRNA", $s1, $e1, qw/. + ./, "ID=${k};Parent=${k}_gene")."\n";
    print {$ofh2} join("\t", $a0->[0], ".", "exon", $s1, $e1, qw/. + ./, "ID=${k}_exon;Parent=$k")."\n";
    print {$ofh2} join("\t", $a0->[0], ".", "CDS",  $s1, $e1, qw/. + ./, "ID=${k}_cds;Parent=$k")."\n";
    $e1 >= $e0 and last; # '=' should be enough.
  }
  $prevEnd = $e0;
}
close ($ofh2);

# &fileSunhh::_rmtree($wrk_dir);

