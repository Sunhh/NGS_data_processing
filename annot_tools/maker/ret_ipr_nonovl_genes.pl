#!/usr/bin/perl -w
# [8/19/2022] Run iprscan to get back good proteins in non-overlap ab init genes.
# [8/23/2022] Fix bugs found in RUNs. Rename output files.
use strict;
use LogInforSunhh;
use fileSunhh;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "help!",
  "ipr_home:s",      # Required: $IPRSCAN_HOME;
  "mkres_dir:s",     # Required: PI596694_SW1;
  "mkpre:s",         # Required: PI596694_SW1.all.maker.non_overlapping_ab_initio.proteins.fasta.gz
  "in_genomeFas:s",  # Required: in.chr.fa
  "opre:s",          # Default: comb; No directories included.
  "pl_retMkGff:s",   # Default: 'perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/maker/ret_maker_abinit_gff3.pl'
  "pl_dealFas:s",    # Default: 'deal_fasta.pl'
  "pl_dealTab:s",    # Default: 'deal_table.pl'
  "pl_dealGff:s",    # Default: 'perl /home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl'
  "exe_mpi:s",       # Default: '/data/Sunhh/install/mpich/bin/mpiexec'
  "binDir_maker:s",  # Default: '/data/Sunhh/src/annotation/maker/maker.3.01.03/bin'
  "mk_cpuN:i",       # Default: 20
  "ipr_cpuN:i",      # Default: 6
  "getFa!",          # Default: Not assigned.
  "noClean!",        # Default: Not assigned. Clean the intermediate files.
);

my $htxt=<<HH;
################################################################################
# perl $0 -ipr_home=\$IPRSCAN_HOME  -mkres_dir PI596694_SW1  -mkpre PI596694_SW1  -in_genomeFas  /path/to/PI596694.chr.fa
#
# -opre          ['comb'] Prefix for iprscan output files.
# -pl_retMkGff   ['perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/maker/ret_maker_abinit_gff3.pl']
# -pl_dealFas    ['deal_fasta.pl']
# -pl_dealTab    ['deal_table.pl']
# -pl_dealGff    ['perl /home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl']
# -exe_mpi       ['/data/Sunhh/install/mpich/bin/mpiexec']
# -binDir_maker  ['/data/Sunhh/src/annotation/maker/maker.3.01.03/bin/']
# -mk_cpuN       [20]
# -ipr_cpuN      [6]
# -getFa         [Boolean] Retrieve CDS and protein fasta files if given.
# -noClean       [Boolean] Do not remove intermediate files if given.
################################################################################
# Example command:
#  perl $0 \\
#   -ipr_cpuN   20 \\
#   -ipr_home   /data/Sunhh/src/Annotation/iprscan/interproscan-5.54-87.0/ \\
#   -mkres_dir  PI596694_SW1/ \\
#   -mkpre      PI596694_SW1 \\
#   -in_genomeFas  ../wmhifi_maker/db/in_genome/PI596694.chr.fa \\
#   -exe_mpi       /data/Sunhh/src/Support/mpich/maker_mpich/install/bin/mpiexec \\
#   -binDir_maker  /data/Sunhh/src/Annotation/maker/maker.3.01.03/bin
#
################################################################################
HH

$opts{'opre'} //= 'comb';
$opts{'pl_retMkGff'}  //= 'perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/maker/ret_maker_abinit_gff3.pl';
$opts{'pl_dealFas'}   //= 'deal_fasta.pl';
$opts{'pl_dealTab'}   //= 'deal_table.pl';
$opts{'pl_dealGff'}   //= 'perl /home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl';
$opts{'exe_mpi'}      //= '/data/Sunhh/install/mpich/bin/mpiexec';
$opts{'binDir_maker'} //= '/data/Sunhh/src/annotation/maker/maker.3.01.03/bin/';
$opts{'mk_cpuN'}      //= 20;
$opts{'ipr_cpuN'}     //= 6;

defined $opts{'help'} and &LogInforSunhh::usage($htxt);

my %inf;
for my $k1 (qw/ipr_home mkres_dir mkpre in_genomeFas/) {
  (defined $opts{$k1} and $opts{$k1} ne "") or &LogInforSunhh::usage("\n-$k1 is required!\n".$htxt);
}
for my $k1 (keys %opts) {
  $inf{$k1} = $opts{$k1};
  unless ($k1 =~ m!^(mkpre|opre|mk_cpuN|ipr_cpuN|help|pl_retMkGff|pl_dealFas|pl_dealTab|pl_dealGff|getFa|noClean)$!i) {
    $inf{$k1} eq '' and die "|$k1|\n";
    $inf{$k1} = &fileSunhh::_abs_path_4link($inf{$k1});
  }
}

# $inf{'in_mkFas'} = "$inf{'mkres_dir'}/$inf{'mkpre'}.all.maker.proteins.fasta";
$inf{'in_mkGff'} = "$inf{'mkres_dir'}/$inf{'mkpre'}_maker.gff3";
$inf{'in_abFas'} = "$inf{'mkres_dir'}/$inf{'mkpre'}.all.maker.non_overlapping_ab_initio.proteins.fasta";
$inf{'in_allGff'} = "$inf{'mkres_dir'}/$inf{'mkpre'}_all.gff3";
$inf{'in_iprDb'}  = "$inf{'ipr_home'}/../ipr.db_slct";
# for my $k1 (qw/in_mkFas in_mkGff in_abFas in_allGff/) {
for my $k1 (qw/in_mkGff in_abFas in_allGff/) {
  -e $inf{$k1} or $inf{$k1} = "$inf{$k1}.gz";
  -e $inf{$k1} or &stopErr("[Err] No maker protein fasta file found for [$k1] [$inf{$k1}]\n");
  $inf{$k1} = &fileSunhh::_abs_path_4link($inf{$k1});
}
$inf{'tDir'} = &fileSunhh::new_tmp_dir('create' => 1);
# $inf{'tDir'} = &fileSunhh::_abs_path_4link($inf{'tDir'});
if ($inf{'in_abFas'} =~ m!\.gz$!) {
  &runCmd("gzip -cd $inf{'in_abFas'} > $inf{'tDir'}/ab.p.fa");
  $inf{'in_abFas'} = &fileSunhh::_abs_path_4link("$inf{'tDir'}/ab.p.fa");
}
$inf{'cur_dir'} = &fileSunhh::_abs_path_4link(".");
$inf{'ipr_dir'} = &fileSunhh::_abs_path_4link("$inf{'mkres_dir'}/iprscan/");
$inf{'genome_pre'} = $inf{'in_genomeFas'}; $inf{'genome_pre'} =~ s!^.+/!!;
$inf{'genome_pre'} =~ s!\.(fa|fas|fasta)$!! or &stopErr("[Err] Bad input fasta file name [$opts{'in_genomeFas'}] [$inf{'genome_pre'}] left\n");

-d $inf{'ipr_dir'} or mkdir($inf{'ipr_dir'}, 0755);
chdir($inf{'ipr_dir'}) or &stopErr("[Err] Failed to get to dir [$inf{'ipr_dir'}]\n");
&runCmd("$inf{'ipr_home'}/interproscan.sh -iprlookup -goterms -pa -etra -cpu $inf{'ipr_cpuN'} -i $inf{'in_abFas'} -b m2 -f tsv");
&runCmd("$inf{'pl_dealTab'} m2.tsv -kSrch_idx $inf{'in_iprDb'} -kSrch_srcCol 3 | cut -f 1 | $inf{'pl_dealTab'} -UniqColLine 0 > m2.tsv.wiIPR.mID");
&runCmd("$inf{'pl_retMkGff'}  m2.tsv.wiIPR.mID  $inf{'in_allGff'} > m2.tsv.wiIPR.0.gff3");
&runCmd("$inf{'binDir_maker'}/maker -CTL");
open F1,'<',"maker_opts.ctl" or die "$!\n";
open O1,'>',"maker_opts.ab.ctl" or die "$!\n";
while (<F1>) {
  chomp;
  if (m!^genome=!) {
    $_ = "genome=$inf{'in_genomeFas'}";
  } elsif (m!^model_org=!) {
    $_ = "model_org=";
  } elsif (m!^repeat_protein=!) {
    $_ = "repeat_protein=";
  } elsif (m!^pred_gff=!) {
    $_ = "pred_gff=m2.tsv.wiIPR.0.gff3";
  } elsif (m!^keep_preds=!) {
    $_ = "keep_preds=1";
  }
  print O1 "$_\n";
}
close O1;
close F1;
&runCmd("$inf{'exe_mpi'} -n $inf{'mk_cpuN'} $inf{'binDir_maker'}/maker maker_opts.ab.ctl maker_bopts.ctl maker_exe.ctl 1> maker.ab-log 2> maker.ab-err");
&runCmd("$inf{'binDir_maker'}/gff3_merge  -d $inf{'genome_pre'}.maker.output/$inf{'genome_pre'}_master_datastore_index.log -g -n -s | grep -v ^# > ab_maker_noFa.gff3");
# &runCmd("$inf{'binDir_maker'}/fasta_merge -d $inf{'genome_pre'}.maker.output/$inf{'genome_pre'}_master_datastore_index.log -o ab");
&runCmd("$inf{'pl_dealGff'} -rmFaInMakerGff -inGff $inf{'in_mkGff'} | grep -v ^# > $inf{'mkres_dir'}/$inf{'opre'}.gff3");
&runCmd("cat ab_maker_noFa.gff3 >> $inf{'mkres_dir'}/$inf{'opre'}.gff3");
if ($opts{'getFa'}) {
  &runCmd("$inf{'pl_dealGff'} -seqret -extractFeat CDS -inGff $inf{'mkres_dir'}/$inf{'opre'}.gff3 -scfFa $inf{'in_genomeFas'} > $inf{'mkres_dir'}/$inf{'opre'}.c.fa");
  &runCmd("$inf{'pl_dealFas'} -cds2aa $inf{'mkres_dir'}/$inf{'opre'}.c.fa | $inf{'pl_dealFas'} -rmDefinition > $inf{'mkres_dir'}/$inf{'opre'}.p.fa");
}
unless ($opts{'noClean'}) {
  &fileSunhh::_rmtree("$inf{'genome_pre'}.maker.output/");
  for my $k2 (qw/maker.ab-err maker.ab-log maker_opts.ab.ctl maker_bopts.ctl maker_evm.ctl maker_exe.ctl maker_opts.ctl temp m2.tsv.wiIPR.0.gff3/) {
    &fileSunhh::_rmtree($k2);
  }
}
chdir($inf{'cur_dir'});
&fileSunhh::_rmtree($inf{'tDir'});

&tsmsg("[Rec] IPR proteins combined for: $inf{'mkpre'}\n");

