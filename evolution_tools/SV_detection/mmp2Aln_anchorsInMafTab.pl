#!/usr/bin/perl
# 5/30/2023 Merge SAM files with small batches.
# 9/8/2023: Add option to control minimap2 parameters.
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;
use fastaSunhh;
use Getopt::Long;
use Parallel::ForkManager;
my %opts;
GetOptions(\%opts,
  "mmp2para:s", # Default: -x asm20 -t 3 -N 20; '-a' is mandatory.
  "help!",
);
$opts{'mmp2para'} //= '-x asm20 -t 3 -N 20';

my $cpuN = 10;

!@ARGV and die "perl $0 out_prefix q.fa r.fa align2w38_anc_fix.maf.blasttab.todo\n";

my $opre = shift;
my $qFas = shift;
my $rFas = shift;
my $mafF = shift;

my $fs_obj = fastaSunhh->new();
my %seqQ = %{$fs_obj->save_seq_to_hash('faFile'=>$qFas)};
my %seqR = %{$fs_obj->save_seq_to_hash('faFile'=>$rFas)};
for (keys %seqQ) {$seqQ{$_}{'seq'} =~ s!\s!!g; $seqQ{$_}{'seq'}=uc($seqQ{$_}{'seq'}); $seqQ{$_}{'len'}=length($seqQ{$_}{'seq'});}
for (keys %seqR) {$seqR{$_}{'seq'} =~ s!\s!!g; $seqR{$_}{'seq'}=uc($seqR{$_}{'seq'}); $seqR{$_}{'len'}=length($seqR{$_}{'seq'});}

my $wd = &fileSunhh::new_tmp_dir('create'=>1);
my @f1 = &fileSunhh::load_tabFile($mafF);

my $pm = new Parallel::ForkManager($cpuN);
for (my $i=0; $i<@f1; $i++) {
  my $pid = $pm->start and next;
  my @ta=@{$f1[$i]};
  my ($qid,$qs,$qe)=@ta[0,6,7];
  my ($rid,$rs,$re)=@ta[1,8,9];
  my $str = '+';
  $qs > $qe and do {$str = '-'; ($qs,$qe)=($qe,$qs); };
  my $qsubSeq = substr($seqQ{$qid}{'seq'}, $qs-1, $qe-$qs+1);
  my $rsubSeq = substr($seqR{$rid}{'seq'}, $rs-1, $re-$rs+1);
  my $qsubID = join("_", $qid, $qs, $qe);
  my $rsubID = join("_", $rid, $rs, $re);
  &fileSunhh::write2file("$wd/q.$i.fa", ">$qsubID\n$qsubSeq\n", '>');
  &fileSunhh::write2file("$wd/r.$i.fa", ">$rsubID\n$rsubSeq\n", '>');
  # &runCmd("minimap2 -a -x asm20 -t 3 -N 20 $wd/r.$i.fa $wd/q.$i.fa | samtools view -h -q 20 > $wd/o.$i.sam");
  &runCmd("minimap2 -a $opts{'mmp2para'} $wd/r.$i.fa $wd/q.$i.fa > $wd/o.$i.sam");
  $pm->finish;
}
$pm->wait_all_children;
my @list_sam_arr = map { "$wd/o.$_.sam" } (0 .. $#f1);
my @list_sam_2;
for (my $i=0; $i*100 < @list_sam_arr; $i++) {
  my $iS = $i*100;
  my $iE = $iS+100-1;
  $iE > $#list_sam_arr and $iE = $#list_sam_arr;
  my $l1 = join(" ", @list_sam_arr[$iS .. $iE]);
  &runCmd("samtools merge -O SAM $wd/m.$i.sam $l1");
  push(@list_sam_2, "$wd/m.$i.sam");
}
if (scalar(@list_sam_2) > 1) {
  my $list_sam_2txt = join(" ", @list_sam_2);
  &runCmd("samtools merge -O SAM ${opre}.sam $list_sam_2txt");
} elsif (scalar(@list_sam_2) == 1) {
  &runCmd("mv $list_sam_2[0] ${opre}.sam");
}

&fileSunhh::_rmtree($wd);

