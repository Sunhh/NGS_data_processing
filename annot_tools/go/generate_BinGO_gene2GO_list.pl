#!/usr/bin/perl
# GO:0003674 MF : molecular_function ; 
# GO:0005575 CC : cellular_component ; 
# GO:0008150 BP : biological_process ; 
use strict; 
use LogInforSunhh; 
use fileSunhh; 

# ==> wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc.sepTbl <==
# eleID	GO_ID	GO_descTxt
# Cla97C01G000010.1	GO:0016021	integral component of membrane
# Cla97C01G000030.1	GO:0004222	metalloendopeptidase activity
# Cla97C01G000030.1	GO:0006508	proteolysis
# 
# ==> ordered_feiID <==
# mrnaID
# Cla97C00G000010.1
# Cla97C00G000020.1
# Cla97C00G000030.1

!@ARGV and die "perl $0 species_ID wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc.sepTbl ordered_feiID > ordered_feiID.gene2GO_BinGO\n"; 

my $ID_species     = shift; 
my $fn_gene2go_sep = shift; 
my $fn_geneList    = shift; 

my %gene2go; 
for ( &fileSunhh::load_tabFile($fn_gene2go_sep, 0) ) {
	my ($geneID, $goID) = @$_; 
	defined $gene2go{$geneID}{'h'}{$goID} and next; 
	$gene2go{$geneID}{'h'}{$goID} = 1; 
	push(@{$gene2go{$geneID}{'a'}}, $goID); 
}
print STDOUT "(species=$ID_species)(type=go)(curator=GO)\n"; 
for ( &fileSunhh::load_tabFile($fn_geneList, 0) ) {
	my $geneID = $_->[0]; 
	$geneID =~ m!^mrnaID$!i and next; 
	for my $topGOID (qw/GO:0003674 GO:0005575 GO:0008150/) {
		defined $gene2go{$geneID}{'h'}{$topGOID} or unshift(@{$gene2go{$geneID}{'a'}}, $topGOID); 
		$gene2go{$geneID}{'h'}{$topGOID} = 1; 
	}
	for my $goID (@{$gene2go{$geneID}{'a'}}) {
		print STDOUT join("\t", $geneID, $goID)."\n"; 
	}
}




