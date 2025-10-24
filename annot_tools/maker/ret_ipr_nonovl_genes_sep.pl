#!/usr/bin/perl -w
# [10/23/2025] Run iprscan to get back good proteins in ab init gene models.
use strict;
use LogInforSunhh;
use fileSunhh;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "help!",
  "ipr_home:s",      # Required: $IPRSCAN_HOME;
  "in_abGff:s",      # Required: Input GFF file for ab initio genes.
  "in_genomeFas:s",  # Required: in.chr.fa
  "opre:s",          # Default: comb; No directories included.
  "pl_dealFas:s",    # Default: 'deal_fasta.pl'
  "pl_dealTab:s",    # Default: 'deal_table.pl'
  "pl_dealGff:s",    # Default: 'perl /home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl'
  "ipr_cpuN:i",      # Default: 6
  "getFa!",          # Default: Not assigned.
  "noClean!",        # Default: Not assigned. Clean the intermediate files.
);

my $htxt=<<HH;
################################################################################
# perl $0 -ipr_home=\$IPRSCAN_HOME  -in_genomeFas  /path/to/PI596694.chr.fa  -in_abGff /path/to/ab_gene.gff3
#
# -opre          ['comb'] Prefix for iprscan output files.
# -pl_dealFas    ['deal_fasta.pl']
# -pl_dealTab    ['deal_table.pl']
# -pl_dealGff    ['perl /home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl']
# -ipr_cpuN      [6]
# -getFa         [Boolean] Retrieve CDS and protein fasta files if given.
# -noClean       [Boolean] Do not remove intermediate files if given.
################################################################################
# Example command:
#  perl $0 \\
#   -ipr_cpuN   20 \\
#   -ipr_home   /data/Sunhh/src/Annotation/iprscan/interproscan-5.54-87.0/ \\
#   -in_genomeFas  ../wmhifi_maker/db/in_genome/PI596694.chr.fa \\
#
################################################################################
HH

$opts{'opre'} //= 'out';
$opts{'pl_dealFas'}   //= 'deal_fasta.pl';
$opts{'pl_dealTab'}   //= 'deal_table.pl';
$opts{'pl_dealGff'}   //= 'perl /home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl';
$opts{'ipr_cpuN'}     //= 6;

defined $opts{'help'} and &LogInforSunhh::usage($htxt);

my %inf;
for my $k1 (qw/ipr_home in_abGff in_genomeFas/) {
  (defined $opts{$k1} and $opts{$k1} ne "") or &LogInforSunhh::usage("\n-$k1 is required!\n".$htxt);
}
for my $k1 (keys %opts) {
  $inf{$k1} = $opts{$k1};
  unless ($k1 =~ m!^(opre|ipr_cpuN|help|pl_dealFas|pl_dealTab|pl_dealGff|getFa|noClean)$!i) {
    $inf{$k1} eq '' and die "|$k1|\n";
    $inf{$k1} = &fileSunhh::_abs_path_4link($inf{$k1});
  }
}
$inf{'in_iprDb'}  = "$inf{'ipr_home'}/../ipr.db_slct";

# Prepare ab initio protein sequence data.
$inf{'tDir'} = &fileSunhh::_abs_path_4link( &fileSunhh::new_tmp_dir('create' => 1) );
$inf{'cur_dir'} = &fileSunhh::_abs_path_4link(".");

$inf{'in_abFas'} = &fileSunhh::_abs_path_4link("$inf{'tDir'}/ab.p.fa");
&runCmd("$inf{'pl_dealGff'} -scfFa $inf{'in_genomeFas'} -inGff $inf{'in_abGff'} -seqret -extractFeat CDS > $inf{'tDir'}/ab.c.fa");
&runCmd("$inf{'pl_dealFas'} -cds2aa $inf{'tDir'}/ab.c.fa|$inf{'pl_dealFas'} -rmTailX_prot|$inf{'pl_dealFas'} -rmDefinition > $inf{'in_abFas'}");

# Identify InterPro domains.
$inf{'ipr_dir'} = &fileSunhh::_abs_path_4link("$inf{'tDir'}/iprscan/");
-d $inf{'ipr_dir'} or mkdir($inf{'ipr_dir'}, 0755);
chdir($inf{'ipr_dir'}) or &stopErr("[Err] Failed to get to dir [$inf{'ipr_dir'}]\n");
&runCmd("$inf{'ipr_home'}/interproscan.sh -iprlookup -goterms -pa -etra -cpu $inf{'ipr_cpuN'} -i $inf{'in_abFas'} -b m2 -f tsv");
&runCmd("$inf{'pl_dealTab'} m2.tsv -kSrch_idx $inf{'in_iprDb'} -kSrch_srcCol 3 | cut -f 1 | $inf{'pl_dealTab'} -UniqColLine 0 > m2.tsv.wiIPR.mID");
&runCmd("$inf{'pl_dealGff'} -inGff $inf{'in_abGff'} -gffret m2.tsv.wiIPR.mID -idType mRNA > m2.tsv.wiIPR.gff3");
chdir($inf{'cur_dir'});

# Generate resulting files.
&fileSunhh::_move( "$inf{'ipr_dir'}/m2.tsv.wiIPR.gff3", "$inf{'opre'}.gff3" );
&fileSunhh::_move( "$inf{'ipr_dir'}/m2.tsv", "$inf{'opre'}.ipr.tsv" );
if ($opts{'getFa'}) {
  &runCmd("$inf{'pl_dealGff'} -seqret -extractFeat CDS -inGff $inf{'opre'}.gff3 -scfFa $inf{'in_genomeFas'} > $inf{'opre'}.c.fa");
  &runCmd("$inf{'pl_dealFas'} -cds2aa $inf{'opre'}.c.fa| $inf{'pl_dealFas'} -rmTailX_prot | $inf{'pl_dealFas'} -rmDefinition > $inf{'opre'}.p.fa");
}
unless ($opts{'noClean'}) {
  &fileSunhh::_rmtree($inf{'tDir'});
}

&tsmsg("[Rec] Proteins with IPR domains kept in: $inf{'opre'}.* files\n");

