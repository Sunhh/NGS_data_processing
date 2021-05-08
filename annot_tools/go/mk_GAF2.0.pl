#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 

# http://www.geneontology.org/page/go-annotation-file-format-20
# The annotation flat file format is comprised of 17 tab-delimited fields.
#
# Column	Content	Required?	Cardinality	Example
# 1	DB	required	1	UniProtKB
# 2	DB Object ID	required	1	P12345
# 3	DB Object Symbol	required	1	PHO3
# 4	Qualifier	optional	0 or greater	NOT
# 5	GO ID	required	1	GO:0003993
# 6	DB:Reference (|DB:Reference)	required	1 or greater	SGD_REF:S000047763|PMID:2676709
# 7	Evidence Code	required	1	IMP
# 8	With (or) From	optional	0 or greater	GO:0000346
# 9	Aspect	required	1	F
# 10	DB Object Name	optional	0 or 1	Toll-like receptor 4
# 11	DB Object Synonym (|Synonym)	optional	0 or greater	hToll|Tollbooth
# 12	DB Object Type	required	1	protein
# 13	Taxon(|taxon)	required	1 or 2	taxon:9606
# 14	Date	required	1	20090118
# 15	Assigned By	required	1	SGD
# 16	Annotation Extension	optional	0 or greater	part_of(CL:0000576)
# 17	Gene Product Form ID	optional	0 or 1	UniProtKB:P12345-2

my $in_goList = shift; 

print STDOUT "!gaf-version: 2.0\n"; 
print STDOUT "!\n!Generated for GO slim\n!\n! http://www.geneontology.org/page/go-annotation-file-format-20\n!\n"; 
my $igo_fh = &openFH($in_goList, '<'); 
while (<$igo_fh>) {
	chomp; 
	my @ta = &splitL("\t", $_); 
	my $class = ( $ta[3] eq "CC" ) ? 'C' : ( $ta[3] eq "MF" ) ? 'F' : ( $ta[3] eq "BP" ) ? 'P' : 'NNNN' ; 
	print STDOUT join("\t", 
		"HS", # database name: UniProtKB
		"$ta[0]", # DB Object ID: P12345
		"$ta[0]", # DB Object Symbol: PHO3
		"", # Qualifier 
		"$ta[1]", # GO ID: GO:0003993
		"HS:0001", # DB:Reference: required, 
		"ISA", # Evidence Code: # 'ND' - No biological Data available; 
		"$ta[1]", # With (or) From
		"$class", # Aspect
		"", # DB Object Name
		"", # DB Object Synonym
		"protein", # DB Object Type
		"taxon:260674", # Taxon 
		"20181129", # Date
		"Swiss-Prot", # Assigned By
		"", # Annotation Extension
		"", # Gene Product Form ID
	)."\n"; 
}
close($igo_fh); 


