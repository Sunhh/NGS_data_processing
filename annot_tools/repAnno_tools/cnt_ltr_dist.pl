#!/usr/bin/perl
use strict;
use warnings;
use fastaSunhh;
use LogInforSunhh;
use fileSunhh;

!@ARGV and die "perl $0 CA02.chr.fa CA02.EDTA.intact.gff3 LTR_dist_est.mao out_prefix\n";

my $exemegacc = 'megacc';
my $exemuscle = 'muscle';
my $plFa2Meg  = 'perl /home/Sunhh/tools/github/NGS_data_processing/reseq_tools/cnvt_tools/fas2meg.pl';

use Parallel::ForkManager;
my $pm = new Parallel::ForkManager(80);

my $faFn  = shift;
my $gffFn = shift;
my $maoFn = shift;
my $opref = shift;

my $fs_obj = fastaSunhh->new();
my %s = %{$fs_obj->save_seq_to_hash('faFile' => $faFn)};
for (keys %s) {
  $s{$_}{'seq'} =~ s![\-\s]!!g;
  $s{$_}{'seq'} = uc($s{$_}{'seq'});
  $s{$_}{'len'} = length($s{$_}{'seq'});
}
my %ltrs;
for my $l1 (&fileSunhh::load_tabFile($gffFn)) {
  $l1->[2] eq 'long_terminal_repeat' or next;
  my $ss = substr($s{$l1->[0]}{'seq'}, $l1->[3]-1, $l1->[4]-$l1->[3]+1);
  $l1->[6] eq '-' and &fastaSunhh::rcSeq(\$ss, 'rc');
  $l1->[8] =~ m!(?:^|;)\s*Parent=(\S+?)\s*(;|$)!i or die "@$l1\n";
  my $pID = $1;
  my $cls = 'NA'; $l1->[8] =~ m!(?:^|;)\s*Classification=(\S+?)\s*(;|$)!i and $cls = $1;
  my $cID = 'NA'; $l1->[8] =~ m!(?:^|;)\s*ID=(\S+?)\s*(;|$)!i and $cID = $1;
  my $ident = -1; $l1->[8] =~ m!(?:^|;)\s*ltr_identity=(\S+?)\s*(;|$)!i and $ident = $1;
  my $motif = 'NA'; $l1->[8] =~ m!(?:^|;)\s*motif=(\S+?)\s*(;|$)!i and $motif = $1;
  my $tsd   = 'NA'; $l1->[8] =~ m!(?:^|;)\s*tsd=(\S+?)\s*(;|$)!i and $tsd = $1;
  my $len   = $l1->[4]-$l1->[3]+1;
  push(@{$ltrs{$pID}}, ["${pID}_${cID}_$l1->[0]_$l1->[3]", $ss, $len, $ident, $cls, $motif, $tsd, "$l1->[0]:$l1->[3]-$l1->[4]:$l1->[6]" ]);
  # ([ltr1: pID_inf, Term_seq, Term_len, Term_ident, class, motif, tsd, Term_pos], [ltr2: ])
}
my @ltrA;
for my $pID (sort keys %ltrs) {
  scalar(@{$ltrs{$pID}}) == 2 or do { &tsmsg("[Wrn] ", scalar(@{$ltrs{$pID}})," terminals in LTR $pID\n"); next;};
  push(@ltrA, [$pID, @{$ltrs{$pID}}]);
}

my $wrkdir = &fileSunhh::new_tmp_dir('create'=>1);
my $absWD  = &fileSunhh::_abs_path($wrkdir);
mkdir("$absWD/out/");
for (my $i=0; $i<@ltrA; $i++) {
  my $pid = $pm->start and next;
  my @a1 = @{$ltrA[$i][1]};
  my @a2 = @{$ltrA[$i][2]};
  &fileSunhh::write2file("$wrkdir/$i.fasta", ">$a1[0]\n$a1[1]\n>$a2[0]\n$a2[1]\n", '>');
  &runCmd("$exemuscle -align $wrkdir/$i.fasta -output $wrkdir/$i.afa");
  &runCmd("$plFa2Meg $wrkdir/$i.afa > $wrkdir/$i.meg");
  &runCmd("$exemegacc -a $maoFn -d $wrkdir/$i.meg -o $absWD/out/$i");
  $pm->finish;
}
$pm->wait_all_children;

open O,'>',"$opref.tab" or die;
print O join("\t", qw/LTR_ID Dist Term1_len Term2_len Term_ident Class Motif TSD Term1_loc Term2_loc/)."\n";
for (my $i=0; $i<@ltrA; $i++) {
  my @a1 = @{$ltrA[$i][1]};
  my @a2 = @{$ltrA[$i][2]};
  open F,'<',"$absWD/out/$i.meg" or die;
  my $dist = 'x';
  while (<F>) {
    m!^\s*\[2\]\s*([eE\-\d.]+)\s*$! and $dist = $1;
  }
  close F;
  print O join("\t", $ltrA[$i][0], $dist, $a1[2], $a2[2], @a1[3..7], $a2[7])."\n";
}
close O;

&fileSunhh::_rmtree($wrkdir);

&tsmsg("[Msg] Done.\n");
