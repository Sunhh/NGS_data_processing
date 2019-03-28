#!/usr/bin/perl
# 2018-11-15 : Edit to use annotated GO terms. 
# GO:0003674 MF : molecular_function ; 
# GO:0005575 CC : cellular_component ; 
# GO:0008150 BP : biological_process ; 
use strict; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"fn_geneList:s", 
	"noAddBasic!", 
	"forBinGO!", 
	"BinGO_speciesID:s", 
); 

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

$opts{'BinGO_speciesID'} //= 'watermelon'; 

my $htxt = <<HH; 
################################################################################
# perl $0 -fn_geneList ordered_feiID  wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc.sepTbl > ordered_feiID.gene2GO_topGO
#
# -noAddBasic         [Boolean] Add CC/BP/MF roots if not given. 
# -forBinGO           [Boolean] Output data for BinGO input. 
#  -BinGO_speciesID   [string] species name to be used. [$opts{'BinGO_speciesID'}]
HH

!@ARGV and &LogInforSunhh::usage($htxt); 
$opts{'help'} and &LogInforSunhh::usage($htxt); 

my $fn_gene2go_sep = shift; 

my %gene2go; 
for ( &fileSunhh::load_tabFile($fn_gene2go_sep, 0) ) {
	my ($geneID, $goID) = @$_; 
	defined $gene2go{$geneID}{'h'}{$goID} and next; 
	$gene2go{$geneID}{'h'}{$goID} = 1; 
	push(@{$gene2go{$geneID}{'a'}}, $goID); 
}
my @geneList; 
if (defined $opts{'fn_geneList'}) {
	@geneList = map { $_->[0] } &fileSunhh::load_tabFile($opts{'fn_geneList'}, 0); 
} else {
	@geneList = sort keys %gene2go; 
}

if ( $opts{'forBinGO'} ) {
	print STDOUT "(species=$opts{'BinGO_speciesID'})(type=go)(curator=GO)\n"; 
}

for my $geneID ( @geneList ) {
	$geneID =~ m!^mrnaID$!i and next; 
	unless ($opts{'noAddBasic'}) {
		for my $topGOID (qw/GO:0003674 GO:0005575 GO:0008150/) {
			defined $gene2go{$geneID}{'h'}{$topGOID} or unshift(@{$gene2go{$geneID}{'a'}}, $topGOID); 
			$gene2go{$geneID}{'h'}{$topGOID} = 1; 
		}
	}
	unless ( defined $gene2go{$geneID}{'h'} ) {
		$gene2go{$geneID}{'a'} = ['NA']; 
	}
	if ( $opts{'forBinGO'} ) {
		for my $goid (@{$gene2go{$geneID}{'a'}}) {
			print join("\t", $geneID, $goid)."\n"; 
		}
	} else {
		print STDOUT join("\t", $geneID, join(",", @{$gene2go{$geneID}{'a'}}))."\n"; 
	}
}




