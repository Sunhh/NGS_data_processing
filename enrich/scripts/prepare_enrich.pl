#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"out:s", 
	"bg_geneList:s", # All genes needed; 
	"in:s", # keggPWYByKO.tab
	"org_colN:i", # 6, column number of organism. 
	"grpID_colN:i", # 2, column number of group ID, default is for Map_ID, and could be for KO_ID. 
	"grpTxt_colN:i", # 3, column number of group text. default is for Map_ID. 
	"dataType:s", # kegg
); 

my %gg; 
&set_Glob(); 

if ( $gg{'dataType'} eq 'kegg' ) {
	&load_background_genelist(); 
	&load_groupAnnot(); 
	&outBG(); 
}


################################################################################
# Sub-routine 
################################################################################
sub set_Glob {
	$gg{'dataType'}     = 'kegg'; 
	$gg{'org_colN'}     = 6; 
	$gg{'grpID_colN'}   = 2; 
	$gg{'grpTxt_colN'}  = 3; 
	$gg{'outFh'}        = \*STDOUT; 
	$gg{'isOfile'}      = 0; 
	$gg{'htxt'} = <<"H1"; 
################################################################################
# perl $0   -bg_geneList genome_geneID_list   -in keggPWYByKO.tab 
#
# -help 
#
#
# -bg_geneList       [filename] All gene IDs in the genome in single-col file. 
# -dataType          ['$gg{'dataType'}'] kegg(go/ipr not included yet)
# -in                [filename] Output of cnvt_keggPWYReconstruction_to_tab.pl ; 
#                      example format : 
#                        PWY_Group1 \\t PWY_Group2               \\t Map_ID \\t Map_description    \\t KO_Num_InMap \\t KO_ID  \\t Org_1 \\n
#                        Metabolism \\t Global and overview maps \\t 01100  \\t Metabolic pathways \\t 831          \\t K00001 \\t Cla97C01G013600.1, Cla97C01G001460.1 \\n
#
#   -org_colN          ['$gg{'org_colN'}'] The organisim column to use. 
#   -grpID_colN        ['$gg{'grpID_colN'}'] column number of group ID, default is for Map_ID, and could be for KO_ID.
#   -grpTxt_colN       ['$gg{'grpTxt_colN'}'] column number of group text. default is for Map_ID. 
#
# -out               [filename] Output to this file if given. 
#
################################################################################
H1
	$opts{'help'} and &LogInforSunhh::usage($gg{'htxt'}); 
	for my $k1 (qw/bg_geneList in/) {
		defined $opts{$k1} or &stopErr("[Err] Lack parameter -$k1\n"); 
		$gg{$k1} = $opts{$k1}; 
	}
	for my $k1 (qw/dataType org_colN grpID_colN grpTxt_colN/) {
		defined $opts{$k1} and $gg{$k1} = $opts{$k1}; 
	}
	$gg{'dataType'} = lc($gg{'dataType'}); 
	if ($gg{'dataType'} eq 'kegg') {
		$gg{'grpID_colN'}  = 2; 
		$gg{'grpTxt_colN'} = 3; 
		for my $k2 (qw/grpID_colN grpTxt_colN/) {
			defined $opts{$k2} and $gg{$k2} = $opts{$k2}; 
		}
	} else {
		&stopErr("[Err] Unknown -dataType $gg{'dataType'}\n"); 
	}
	defined $opts{'out'} and do { $gg{'outFh'} = &openFH($opts{'out'}, '>'); $gg{'isOfile'} = 1; }; 
	return; 
}# set_Glob() 

sub load_background_genelist {
	# Load file -bg_geneList into 
	#  $gg{'data'}{'bg_eleID_h'} = { eleID->order_number }; 
	#  $gg{'data'}{'bg_eleID_a'} = [ eleID, eleID, ... ]; # non-redundant 
	#
	my $fh = &openFH($gg{'bg_geneList'}); 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		defined $gg{'data'}{'bg_eleID_h'}{$ta[0]} and next; 
		$gg{'data'}{'bg_eleID_h'}{$ta[0]} = $.; 
		push(@{$gg{'data'}{'bg_eleID_a'}}, $ta[0]); 
	}
	close($fh); 
	return; 
}# load_background_genelist() 

sub load_groupAnnot {
	# Load -in file into 
	#   $gg{'data'}{'grpID2Txt'} = { grpID=>grpAnnotTxt }
	#   $gg{'data'}{'ele2grp'}   = { eleID=>[grpID, grpID, ...] }
	my $fh = &openFH($gg{'in'}, '<'); 
	my %has_ele2grp; 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		$ta[$gg{'org_colN'}] =~ m!^(\s*|NA)$!i and next; 
		my $grpID = $ta[$gg{'grpID_colN'}]; 
		my $grpTxt = $ta[$gg{'grpTxt_colN'}]; 
		$gg{'data'}{'grpID2Txt'}{$grpID} //= $grpTxt; 
		$gg{'data'}{'grpID2Txt'}{$grpID} eq $grpTxt or &stopErr("[Err] Different annotation for [$grpID]: '$gg{'data'}{'grpID2Txt'}{$grpID}' VS. '$grpTxt'\n"); 
		my @gen = map { s!^\s*|\s*$!!g; $_; } &splitL(",", $ta[$gg{'org_colN'}]); 
		for my $g1 (@gen) {
			defined $gg{'data'}{'bg_eleID_h'}{$g1} or do { &tsmsg("[Wrn] Ignore unknown element ID $g1 in $ta[$gg{'org_colN'}]\n"); next; }; 
			defined $has_ele2grp{$g1}{$grpID} and next; 
			$has_ele2grp{$g1}{$grpID} = 1; 
			push(@{$gg{'data'}{'ele2grp'}{$g1}}, $grpID); 
		}
	}
	close($fh); 
	return; 
}# load_groupInfo() 


sub outBG {
	print {$gg{'outFh'}} join("\t", qw/eleID grpID grpDescription/)."\n"; 
	for my $g1 (@{$gg{'data'}{'bg_eleID_a'}}) {
		if (defined $gg{'data'}{'ele2grp'}{$g1}) {
			for my $a1 (@{$gg{'data'}{'ele2grp'}{$g1}}) {
				print {$gg{'outFh'}} join("\t", $g1, $a1, $gg{'data'}{'grpID2Txt'}{$a1})."\n"; 
			}
		} else {
			print {$gg{'outFh'}} join("\t", $g1, 'NA', 'NA')."\n"; 
		}
	}
	$gg{'isOfile'} and close($gg{'outFh'}); 
}# outBG() 

