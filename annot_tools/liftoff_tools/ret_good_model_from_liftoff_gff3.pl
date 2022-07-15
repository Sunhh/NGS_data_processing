#!/usr/bin/perl
# 5/24/2022: Filter out bad gene models mapped by Liftoff.
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;

my $pl_deal_gff3 = 'perl /home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl';
!@ARGV and die "perl $0 liftoff_out.gff3 > filtered.gff3\n";

my $inFn = shift;

my $wdir = &fileSunhh::new_tmp_dir('create' => 1);

open F,'<',"$inFn" or die;
open OM,'>',"$wdir/good_mrna" or die;
open OG,'>',"$wdir/good_gene" or die;

while (<F>) {
  m!^\s*(#|$)! and next;
  my @ta=split(/\t/, $_);
  if ($ta[2] =~ m!^gene$!i) {
    $ta[8] =~ m!(low_identity=True|valid_ORFs=0|partial_mapping=True)! and next;
    $ta[8] =~ m!^ID=([^\s;]+)! or &stopErr("[Err] Failed to parse ID from liftoff gff3 in lin: $_\n");
    print OG "$1\n";
  } elsif ($ta[2] =~ m!^mRNA$!i) {
    $ta[8] =~ m!(partial_ORF=True|inframe_stop_codon=True|missing_stop_codon=True|missing_start_codon=True|valid_ORF=False)! and next;
    $ta[8] =~ m!^ID=([^\s;]+)! or &stopErr("[Err] Failed to parse ID from liftoff gff3 in lin: $_\n");
    print OM "$1\n";
  }
}
close OM;
close OG;
&runCmd("$pl_deal_gff3 -inGff $inFn -gffret $wdir/good_gene -idType 'gene' > $wdir/good_gene.gff3");
&runCmd("$pl_deal_gff3 -inGff $wdir/good_gene.gff3 -gffret $wdir/good_mrna -idType 'mRNA' > $wdir/good_mrna.gff3");
open F2,'<',"$wdir/good_mrna.gff3" or die ;
while (<F2>) {
  print STDOUT $_;
}
close F2;

&fileSunhh::_rmtree($wdir);

