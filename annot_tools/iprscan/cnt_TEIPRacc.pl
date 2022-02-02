#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;

!@ARGV and die "perl $0 potential_TE_IPRacc ipr_all6_tsv.TEprot.IPRacc.line > ipr_all6_tsv.TEprot.IPRacc.line.cnt\n";

my $f1 = shift;
my $f2 = shift;

my @a1 = &fileSunhh::load_tabFile($f1);
my @a2 = &fileSunhh::load_tabFile($f2);

my %teID; 
for my $l1 (@a1) {
  $teID{$l1->[0]} = 1;
}
print join("\t", qw/geneID Is_TE ttl_IPRacc TE_IPRacc TE_IPRacc_Perc IPRaccessions IPRannotations/)."\n";
for my $l2 (@a2) {
  my @iprIDs = split(/;;/, $l2->[1]);
  my $cnt_ttl = scalar(@iprIDs);
  my $cnt_te  = 0;
  for my $id (@iprIDs) {
    defined $teID{$id} and $cnt_te ++;
  }
  my $is_te = 'notTE';
  $cnt_ttl == $cnt_te and $is_te = 'TE';
  print join("\t", $l2->[0], $is_te, $cnt_ttl, $cnt_te, sprintf("%.2f", 100*$cnt_te/$cnt_ttl), @{$l2}[1..$#$l2])."\n";
  
}

# ==> potential_TE_IPRacc <==
# IPR012337	Ribonuclease H-like superfamily
# IPR043502	DNA/RNA polymerase superfamily
# IPR036397	Ribonuclease H superfamily
# IPR021109	Aspartic peptidase domain superfamily
# IPR043128	Reverse transcriptase/Diguanylate cyclase domain
# IPR000477	Reverse transcriptase domain
# IPR002156	Ribonuclease H domain
# IPR013103	Reverse transcriptase, RNA-dependent DNA polymerase
# IPR005162	Retrotransposon gag domain
# IPR018289	MULE transposase domain
# 
# ==> ipr_all6_tsv.TEprot.IPRacc.line <==
# CaU00G01120.1	IPR016024	Armadillo-type fold
# CaU00G01180.1	IPR012337;;IPR023211;;IPR036397;;IPR043502;;IPR044925	Ribonuclease H-like superfamily;;DNA polymerase, palm domain superfamily;;Ribonuclease H superfamily;;DNA/RNA polymerase superfamily;;His-Me finger superfamily
# CaU00G01210.1	IPR012337	Ribonuclease H-like superfamily
# CaU00G01420.1	IPR005135;;IPR036691	Endonuclease/exonuclease/phosphatase;;Endonuclease/exonuclease/phosphatase superfamily
# CaU00G01430.1	IPR000850;;IPR006266;;IPR027417;;IPR033690	Adenylate kinase/UMP-CMP kinase;;UMP-CMP kinase;;P-loop containing nucleoside triphosphate hydrolase;;Adenylate kinase, conserved site
# CaU00G01460.1	IPR001951;;IPR004823;;IPR005135;;IPR009072;;IPR019809;;IPR035425;;IPR036691	Histone H4;;TATA box binding protein associated factor (TAF);;Endonuclease/exonuclease/phosphatase;;Histone-fold;;Histone H4, conserved site;;CENP-T/Histone H4, histone fold;;Endonuclease/exonuclease/phosphatase superfamily
# CaU00G01520.1	IPR036691	Endonuclease/exonuclease/phosphatase superfamily
# CaU00G01540.1	IPR001965;;IPR004875;;IPR011011;;IPR013083;;IPR019787;;IPR036397	Zinc finger, PHD-type;;DDE superfamily endonuclease domain;;Zinc finger, FYVE/PHD-type;;Zinc finger, RING/FYVE/PHD-type;;Zinc finger, PHD-finger;;Ribonuclease H superfamily
# CaU00G01590.1	IPR001584;;IPR012337;;IPR036397;;IPR036444	Integrase, catalytic core;;Ribonuclease H-like superfamily;;Ribonuclease H superfamily;;Phospholipase A2 domain superfamily
# CaU00G01660.1	IPR001878	Zinc finger, CCHC-type
