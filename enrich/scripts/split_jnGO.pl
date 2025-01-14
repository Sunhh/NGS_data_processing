#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"in:s", 
	"col_geneID:i", 
	"col_jnGO:i", 
	"skipHlineN:i", 
	"help!", 
); 

$opts{'col_geneID'} //= 0; 
$opts{'col_jnGO'}   //= 2; 
$opts{'skipHlineN'} //= 0; 

my $htxt = <<"H1"; 
################################################################################
# perl $0 -in wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc > wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc.2col
#
# -help
#
# -col_geneID       [$opts{'col_geneID'}] Column-N of geneID ; 
# -col_jnGO         [$opts{'col_jnGO'}] Column-N of joined GO_txt; 
#                      GO_txt example : 
#                      GO:0004222 (metalloendopeptidase activity), GO:0006508 (proteolysis), GO:0016021 (integral component of membrane)
#
# -skipHlineN       [$opts{'skipHlineN'}] Number of lines to be skipped as header. 
#
#
# OUT example : 
#   Cla97C01G000030.1 \\t GO:0004222 \\t metalloendopeptidase activity \\n
#   Cla97C01G000030.1 \\t GO:0006508 \\t proteolysis \\n
#   Cla97C01G000030.1 \\t GO:0016021 \\t integral component of membrane \\n
#
################################################################################
# Example of input file wm97pbV2ID_feiID_clean.prot.b2g_wiIPR.annot.all_desc : 
# mrnaID	GO_funcDesc	GO_IDdesc	GO_EC
# Cla97C01G000010.1	mitochondrial import inner membrane translocase subunit tim22-like	GO:0016021 (integral component of membrane)	NA
# Cla97C01G000030.1	membrane-bound transcription factor site-2 protease homolog isoform x1	GO:0004222 (metalloendopeptidase activity), GO:0006508 (proteolysis), GO:0016021 (integral component of membrane)	EC:3.4.24.0
# Cla97C01G000040.1	npl4-like protein 1	GO:0005634 (nucleus), GO:0006511 (ubiquitin-dependent protein catabolic process), GO:0031625 (ubiquitin protein ligase binding), GO:0043130 (ubiquitin binding)	NA
# 
# Example of result file : 
# Cla97C01G000010.1	GO:0016021	integral component of membrane
# Cla97C01G000030.1	GO:0004222	metalloendopeptidase activity
# Cla97C01G000030.1	GO:0006508	proteolysis
# Cla97C01G000030.1	GO:0016021	integral component of membrane
# Cla97C01G000040.1	GO:0005634	nucleus
# Cla97C01G000040.1	GO:0006511	ubiquitin-dependent protein catabolic process
# Cla97C01G000040.1	GO:0031625	ubiquitin protein ligase binding
# Cla97C01G000040.1	GO:0043130	ubiquitin binding
# Cla97C01G000050.1	GO:0016787	hydrolase activity
# Cla97C01G000060.1	GO:0016021	integral component of membrane
################################################################################
H1

$opts{'help'} and &LogInforSunhh::usage($htxt); 
defined $opts{'in'} or &LogInforSunhh::usage($htxt); 

my $fh = &openFH("$opts{'in'}", '<'); 
print STDOUT join("\t", qw/eleID GO_ID GO_descTxt/)."\n"; 
JNGOLINE:
while (<$fh>) {
	$. <= $opts{'skipHlineN'} and next; 
	chomp; 
	my @ta = &splitL("\t", $_); 
	my $geneID = $ta[$opts{'col_geneID'}]; 
	my $jnGO_txt = $ta[$opts{'col_jnGO'}]; 
	my @sepGO_txt; 
	while ($jnGO_txt =~ s!^\s*,*\s*(GO:\d+)\s*!!) {
		push(@sepGO_txt, [$1]); 
		if ($jnGO_txt =~ m!^,!) {
			$jnGO_txt =~ s!^,\s*!!; 
			$sepGO_txt[-1][1] = 'NA'; 
		} elsif ($jnGO_txt =~ s!^\(!!) {
			my $followTxt = ''; 
			my $i_brack = 1; 
			for (my $i=0; $i<length($jnGO_txt); $i++) {
				my $curChar = substr($jnGO_txt, $i, 1); 
				if ($curChar eq '(') {
					$i_brack ++; 
				} elsif ($curChar eq ')') {
					$i_brack --; 
				}
				if ($i_brack > 0) {
					$followTxt .= $curChar; 
				} else {
					substr($jnGO_txt, 0, $i+1) = ''; 
					last; 
				}
			}
			$sepGO_txt[-1][1] = $followTxt; 
		} else {
			&tsmsg("[Err][Wrn] Please manually parse jnGO_txt section for line : [$_] Resting [$jnGO_txt]\n"); 
			print STDOUT join("\t", $geneID, $ta[$opts{'col_jnGO'}])."\n"; 
			next JNGOLINE; 
		}
	}
	unless ( $jnGO_txt =~ m!^\s*$! ) {
		&tsmsg("[Err][Wrn] Please manually parse jnGO_txt section for line : [$_] Resting [$jnGO_txt]\n"); 
		print STDOUT join("\t", $geneID, $ta[$opts{'col_jnGO'}])."\n"; 
		next JNGOLINE; 
	}
	for my $t1 (@sepGO_txt) {
		print STDOUT join("\t", $geneID, @{$t1}[0,1])."\n"; 
	}
}
close($fh); 


