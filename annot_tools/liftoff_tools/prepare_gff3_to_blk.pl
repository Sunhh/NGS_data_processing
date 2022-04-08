#!/usr/bin/perl
use strict;
use warnings;
use LogInforSunhh;

for my $tag (qw/CL CM CA CC/) {
  &runCmd("perl /home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl -inGff input/${tag}pan.trim2CDS.gff3 -getJnLoc > input/${tag}pan.trim2CDS.gff3.JnLoc");
  open F,'<',"input/${tag}pan.trim2CDS.gff3.JnLoc" or die;
  open O,'>',"input/${tag}pan.trim2CDS.blk" or die;
  while (<F>) {
    chomp;
    my @ta=split(/\t/, $_);
    $ta[0] eq 'mrnaID' and next;
    print O join("\t", "$tag:$ta[0]", @ta[2,5,9])."\n";
  }
  close O;
  close F;
  &runCmd("rm input/${tag}pan.trim2CDS.gff3.JnLoc");
}

