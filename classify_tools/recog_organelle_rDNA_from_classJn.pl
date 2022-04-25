#!/usr/bin/perl
# [4/9/2022] Need to convert this as pipeline some day later for repeat works!!!
# [4/15/2022] Merge chloroplast and plastid to chloroplast class.
# [4/15/2022] Try to recognize bacteria classes.
# [4/15/2022] Try to identify dispersed hits for rDNAs and satellites.
# [4/21/2022] Make this specific for Euk as In class.
# [4/22/2022] Allow more chloroplast
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;

!@ARGV and die "perl $0 HAl2u.noRed.tochk.toNt.bn6.highCopy.class.jn > HAl2u.noRed.tochk.toNt.bn6.highCopy.class.jn.s1\n";

my $minPerc_organelle = 0.95;
my $minPerc_Ex        = 0.50;
my $minPerc_strict    = 0.95;
my $minPerc_microorg  = 0.50;
my $minPerc_TE        = 0.50;
my $maxInLen_type1    = 1000;
my $maxInPerc_type1   = 0.05;

while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[1] eq 'SeqLen' and next;
  # To identify rDNA and satellites:
  if      ( $ta[5] =~ m!^rDNA:\d+(;;(Eukaryota|Satellite):\d+)*$!i ) {
    if      ( $ta[4]+$ta[8] >= $ta[1] * $minPerc_TE and ($ta[7] < $maxInLen_type1 or $ta[7] < $ta[1] * $maxInPerc_type1) ) {
      # I'd like to use a smaller ExMatLen/SeqLen ratio to identify rDNA related region.
      # I'd like to require no more than 1 kb hits to In class to assure there are no genes in this contig.
      print join("\t", "rDNA", $_)."\n"; next;
    }
  } elsif ( $ta[5] =~ m!^Satellite:\d+(;;(Eukaryota|rDNA):\d+)*$!i ) {
    if ( $ta[4]+$ta[8] >= $ta[1] * $minPerc_TE and ($ta[7] < $maxInLen_type1 or $ta[7] < $ta[1] * $maxInPerc_type1) ) {
      print join("\t", "Satellite", $_)."\n"; next;
    }
  }
  # To identify bacteria and viruses.
  ### I'd like to use a similar way as rDNA and satellites.
  if      ($ta[5] =~ m!^Viruses:\d+$!i) {
    if ( $ta[4] >= $ta[1] * $minPerc_microorg and $ta[7] < $maxInLen_type1 ) {
      print join("\t", "Viruses", $_)."\n"; next;
    }
  } elsif ($ta[5] =~ m!^Bacteria:\d+$!i) {
    if ( $ta[4] >= $ta[1] * $minPerc_microorg and $ta[7] < $maxInLen_type1 ) {
      print join("\t", "Bacteria", $_)."\n"; next;
    }
  } elsif ($ta[5] =~ m!^Archaea:\d+$!i) {
    if ( $ta[4] >= $ta[1] * $minPerc_microorg and $ta[7] < $maxInLen_type1 ) {
      print join("\t", "Archaea", $_)."\n"; next;
    }
  } elsif ($ta[5] =~ m!^NA:\d+$!i) {
    if ( $ta[4] >= $ta[1] * $minPerc_microorg and $ta[7] < $maxInLen_type1 ) {
      print join("\t", "NA", $_)."\n"; next;
    }
  } elsif ($ta[5] =~ m!^((?:NA|Viruses|Bacteria|Archaea):\d+(?:;;)?)+$!i) {
    if ( $ta[4] >= $ta[1] * $minPerc_microorg and $ta[7] < $maxInLen_type1 ) {
      print join("\t", "NA", $_)."\n"; next;
    }
  }

  # To identify organelle genome and others.
  $ta[4] >= $ta[1] * $minPerc_Ex or next; # This is loose to allow sharing between nuclear genome and mtDNA/ctDNA.
  if      ($ta[5] =~ m!^((?:Chloroplast|Plastid):\d+(?:;;)?)+(;;Eukaryota:\d+)?$!i and ($ta[4]+$ta[8]) >= $ta[1] * $minPerc_organelle) {
    # Allow shared region between nuclear and mt/ct genomes.
    print join("\t", "Chloroplast", $_)."\n";
  } elsif ($ta[5] =~ m!^Mitochondrion:\d+(;;Eukaryota:\d+)?$!i and ($ta[4]+$ta[8]) >= $ta[1] * $minPerc_organelle) {
    # Allow shared region between nuclear and mt/ct genomes.
    print join("\t", "Mitochondrion", $_)."\n";
  } elsif ($ta[5] =~ m!^rDNA:\d+(;;Eukaryota:\d+)?$!i and ($ta[4]+$ta[8]) >= $ta[1] * $minPerc_strict) {
    # Allow shared region between nuclear and mt/ct genomes.
    print join("\t", "rDNA", $_)."\n";
  } elsif ($ta[5] =~ m!^Satellite:\d+(;;Eukaryota:\d+)?$!i and ($ta[4]+$ta[8]) >= $ta[1] * $minPerc_strict) {
    # Allow shared region between nuclear and mt/ct genomes.
    print join("\t", "Satellite", $_)."\n";
  } elsif ($ta[4] >= $ta[1] * $minPerc_strict) {
    if      ($ta[5] =~ m!^Viruses:\d+$!i) {
      print join("\t", "Viruses", $_)."\n";
    } elsif ($ta[5] =~ m!^Bacteria:\d+$!i) {
      print join("\t", "Bacteria", $_)."\n";
    } elsif ($ta[5] =~ m!^Archaea:\d+$!i) {
      print join("\t", "Archaea", $_)."\n";
    } elsif ($ta[5] =~ m!^NA:\d+$!i) {
      print join("\t", "NA", $_)."\n";
    } elsif ($ta[5] =~ m!^((?:NA|Viruses|Bacteria|Archaea):\d+(?:;;)?)+$!i) {
      print join("\t", "NA", $_)."\n";
    } else {
      # This ratio can be stricter or looser. To me 95% coverage is safe enough to lower the false positive.
      my @b = sort { $b->[1] <=> $a->[1] } map { [split(/:/, $_)] } split(/;;/, $ta[5]);
      if      ($b[0][0] =~ m!^(Chloroplast|Mitochondrion|Plastid)$!i) {
        print join("\t", "Organelle", $_)."\n";
      } elsif ($b[0][0] =~ m!^Satellite$!i) {
        print join("\t", "Satellite", $_)."\n";
      } elsif ($b[0][0] =~ m!^rDNA$!i) {
        print join("\t", "rDNA", $_)."\n";
      } else {
        # Should not reach here!
        print join("\t", "Others", $_)."\n";
      }
    }
  } else {
    # This is not assigned as contamination.
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


