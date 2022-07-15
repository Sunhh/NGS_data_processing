#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;
use LogInforSunhh;

!@ARGV and die "perl $0  final_C39.ahrd  C31.ndf_snps_rmGap.ref.effG.slct1 > C31.ndf_snps_rmGap.ref.effG.slct1_ahrd\n";

my $fn_ahrd = shift;

my %g2ahrd = map { $_->[0] => $_->[1] } &fileSunhh::load_tabFile($fn_ahrd, '<');

while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'VAR_ID' and do { print "$_\n"; next; };
  $ta[5] eq '-' and do { print "$_\n"; next; };
  my @genes;
  for my $tb (split(/;/, $ta[5])) {
    defined $g2ahrd{$tb} or &stopErr("[Err] failed to find AHRD for gene [$tb]\n");
    push(@genes, "$tb - $g2ahrd{$tb}");
  }
  $ta[5] = join(";", @genes);
  print join("\t", @ta)."\n";
}

# ==> C31.ndf_snps_rmGap.qry.effG.slct1 <==
# VAR_ID  Chr_ID  Chr_start       Chr_end Affected_type   Affected_genes  VAR_type1       VAR_type2       VAR_len VAR_annot
# SNP_300 C39_Chr01       6775    6794    -       -       DELETION        deletion        20      ID=SNP_300;Name=deletion;del_len=20;query_dir=1;query_sequence=C31_Chr01;query_coord=16263;query_bases=-;ref_bases=GGATCCCGCTTGGCTGCTGC;color=#0000EE
# SNP_303 C39_Chr01       6965    6965    -       -       INSERTION       insertion       235     ID=SNP_303;Name=insertion;ins_len=235;query_dir=1;query_sequence=C31_Chr01;query_coord=16435-16669;query_bases=GGTCAATTCAGGGTTTCTTTTTCTTCCGAGGAAAAGCGATGGGCGGCGGTCGGGGGACAGAATGAGCTAGAGATAAGGAATTTATATGCAAGGGGACCATGAAATCGAACTTGGGACATTGTAGAAATACCTCCAGGACTAAAATCAAGTTTCATAGTACTCCTTTCTAGGAGCTGAATCCTCCAACGTAAGAACTTGGTCCGCTTTTCCTAATACTCCTCTCCCTGCTGTGTCC;ref_bases=-;color=#EE0000
# SNP_317 C39_Chr01       8432    8432    -       -       INSERTION       insertion       698     ID=SNP_317;Name=insertion;ins_len=698;query_dir=1;query_sequence=C31_Chr01;query_coord=18128-18825;query_bases=AGATAGAGATAGTTCAAAGAACTTGCTTTCTTTCTCCTTCTCTTCTTCTTCTCTTCGATGTCCCTCTTGTTTAACAGTAGCAGATTCGTTCACTGCAACGCGGTCTGATTCAGAATTAGGCATCATCTTTATCTTTGGATGCTCCATTCCCTCTCTATCTCCATGGCAATCATAAGAATAGTATCGATTCAAACTGTTCCTATATATTATTTGTTCCTCGCCCCAGCTCCAGATTTTTGCGCTGATATGATCAGAATATGAGAAGAAAACATGTTAATAGAGACAAAGAAGGGAAAAACAATGACAAATGAATACCCAATGTTTGTTGTTACTTTAGGCACCTTTTTGACGCTGAAATTTGAGCTGTTCTTGATGTATACATTATTCTGCCTCCCTTTGCTACTTCGCCCAAATCCCTCACGTTTTCTTCTTCTACGCACCAGTCGTCTTCGTCTTCCTCCTCGTCCTCGAAGATACTACGCTTTCTTTTAAGGGTTGAAGAATTTAGATCTGTGAACCGAACAAAAAATCATTATTGAAAAGAAAACAAAATTTACGCCGGAAATGTGTGTACAGGAAAATTGCTAATATTAGGAGAGTTGAAAGGGCAATAGCTTTTTGTGCCATTTACAACAATTTGTGTATGAACTAAAAATAATAAAACATTGGGTTAAATTTGAAATAATAAAGAGGGTTAT;ref_bases=-;color=#EE0000
# SNP_348 C39_Chr01       11271   11271   intron  C39C01G000010.1 INSERTION       insertion       24      ID=SNP_348;Name=insertion;ins_len=24;query_dir=1;query_sequence=C31_Chr01;query_coord=21662-21685;query_bases=CATTTATTGGGGAAGAAATTCAGC;ref_bases=-;color=#EE0000
# SNP_381 C39_Chr01       12981   13030   intron  C39C01G000010.1 DELETION        deletion        50      ID=SNP_381;Name=deletion;del_len=50;query_dir=1;query_sequence=C31_Chr01;query_coord=23419;query_bases=-;ref_bases=TTAGTACCGTTTTGACACCCAATAAGAGCTTATGCTTCCAGAGCCCAAGC;color=#0000EE
# SNP_383 C39_Chr01       13677   13677   CDS     C39C01G000010.1 INSERTION       insertion       20      ID=SNP_383;Name=insertion;ins_len=20;query_dir=1;query_sequence=C31_Chr01;query_coord=24069-24088;query_bases=ACACTTATTAAGAGGGATCA;ref_bases=-;color=#EE0000
# SNP_395 C39_Chr01       14700   14733   upstream        C39C01G000020.1;C39C01G000010.1 DELETION        deletion        34      ID=SNP_395;Name=deletion;del_len=34;query_dir=1;query_sequence=C31_Chr01;query_coord=25133;query_bases=-;ref_bases=CTCTCTCTCTCTCTCTTTCTCTCTCTCTCTCTCC;color=#0000EE
# SNP_446 C39_Chr01       21521   21544   -       -       DELETION        deletion        24      ID=SNP_446;Name=deletion;del_len=24;query_dir=1;query_sequence=C31_Chr01;query_coord=31936;query_bases=-;ref_bases=TTTTTTTTTTTTTTTTTATATTTA;color=#0000EE
# SNP_448 C39_Chr01       21584   21584   -       -       INSERTION       insertion       373     ID=SNP_448;Name=insertion;ins_len=373;query_dir=1;query_sequence=C31_Chr01;query_coord=31979-32351;query_bases=TAAGCTTGTATGTTTGCTCATTCCTAACATACAAGAGATGCGTTGTTCAATATCTAGAATACTACCAAAGTCATTTTGCTTACTGCTTTTTTCTTCTTGTTTGGAAAGATAGCATGCTAGTTTTTAGACTTGATAGAAGGAAAAAAATGATCTTTTTATTTATTTCTGTAGGCTTGTAAAAAAATTCTTCAACCAATGCTAAAGGGCGATTGGAAAAAAAAAATGCAATAGAATGTAATAATTATTACGATATAATACTTATTTTTTACTAATTATATTTATTATTAGTATTTTTAATTTTAATTTCAATCAATTAATTATATTTTAAAATAAAATACATTTAACTAAAAGAAAAAGAAAAAGAAAAAAAAGG;ref_bases=-;color=#EE0000

