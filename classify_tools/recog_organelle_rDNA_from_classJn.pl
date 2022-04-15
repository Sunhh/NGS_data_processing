#!/usr/bin/perl
# [4/9/2022] Need to convert this as pipeline some day later for repeat works!!!
# [4/15/2022] Merge chloroplast and plastid to chloroplast class.
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;

!@ARGV and die "perl $0 HAl2u.noRed.tochk.toNt.bn6.highCopy.class.jn > HAl2u.noRed.tochk.toNt.bn6.highCopy.class.jn.s1\n";

while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[1] eq 'SeqLen' and next;
  $ta[4] >= $ta[1] * 0.95 or next; # This ratio can be stricter or looser. To me 95% coverage is safe enough to lower the false positive.
  if      ($ta[5] =~ m!^rDNA:\d+$!i) {
    print join("\t", "rDNA", $_)."\n";
  } elsif ($ta[5] =~ m!^((?:Chloroplast|Plastid):\d+(?:;;)?)+$!i) {
    print join("\t", "Chloroplast", $_)."\n";
  } elsif ($ta[5] =~ m!^Mitochondrion:\d+$!i) {
    print join("\t", "Mitochondrion", $_)."\n";
  } elsif ($ta[5] =~ m!^Satellite:\d+$!i) {
    print join("\t", "Satellite", $_)."\n";
  } else {
    my @b = sort { $b->[1] <=> $a->[1] } map { [split(/:/, $_)] } split(/;;/, $ta[5]);
    if      ($b[0][0] =~ m!^(Chloroplast|Mitochondrion|Plastid)$!i) {
      print join("\t", "Organelle", $_)."\n";
    } elsif ($b[0][0] =~ m!^Satellite$!i) {
      print join("\t", "Satellite", $_)."\n";
    } else {
      # Should not reach here!
      print join("\t", "Others", $_)."\n";
    }
  }
}


# deal_fasta.pl -keep_len 0-1000000 ../db/C31.hf2a.noRed.fa > C31.tochk.fa
# blastn -query C31.tochk.fa -out C31.tochk.fa.toNt.bn6 \
#  -db nt -evalue 1e-5 -num_threads 50 -max_hsps 50 -max_target_seqs 50 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand staxids sscinames sskingdoms stitle'
# perl /home/Sunhh/tools/github/NGS_data_processing/classify_tools/classify_region_byBn6.pl \
#  C31.tochk.fa.toNt.bn6 \
#  -joinInEx   C31.tochk.fa.toNt.bn6.jnInEx \
#  -InList   Eukaryota:Satellite \
#  -ExList   NA:Viruses:Bacteria:Archaea:rDNA:Chloroplast:Mitochondrion:Plastid \
# > C31.tochk.fa.toNt.bn6.class
# perl /home/Sunhh/tools/github/NGS_data_processing/classify_tools/get_Ex_region.pl C31.tochk.fa.toNt.bn6.class 1> C31.tochk.fa.toNt.bn6.class.sep 2> C31.tochk.fa.toNt.bn6.class.jn
# awk 'NR != 1 { print $5/$2*100"\t"$0 }' C31.tochk.fa.toNt.bn6.class.jn > C31.tochk.fa.toNt.bn6.class.jn.1
# cat C31.tochk.fa.toNt.bn6.class.jn.1 | grep rDNA  | awk '{print $1"\t"$2"\trDNA"}'    >  C31.tochk.fa.toNt.bn6.class.jn.2
# cat C31.tochk.fa.toNt.bn6.class.jn.1 | grep -v rDNA | awk '{print $1"\t"$2"\torganelle"}' >> C31.tochk.fa.toNt.bn6.class.jn.2


