#!/usr/bin/perl
use strict; 
use warnings; 
use Time::Piece; 
use fileSunhh; 
use LogInforSunhh; 
use XML::LibXML; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"task:s", # convert 
	  "dbXml2GOlist!", "xpath2go:s", 
	  "xml5_to_raw4!", 
	  "join_xml5!", 
	  "dbV5_to_dbV4!", 
	   "dirV5:s", "dirV4old:s", "dirV4new:s", 
	   "hmmconvert:s", "doIndex!", 
	   "f_matXML:s", 
	   "f_matDTD:s", 
	   "f_iprXML:s", 
	   "f_iprDTD:s", 
	   "f_superF:s", 
	   "f_smrtTh:s", 
	   "f_smrtDe:s", 
	   "f_smrtHm:s", 
	   "f_sfHMM:s", 
	   "f_prints:s", 
	   "f_gene3d:s",    
	   "f_tigrfam:s", 
	   "f_pfamAH:s", 
	   "f_pfamAS:s", 
	   "f_pfamC:s",     
	   "f_prodom:s", 
	   "f_panther:s", 
	   "f_prosEv:s", 
	   "f_prosite:s", 
	   "f_fprint:s", 
	   "f_superFT:s", 
	   "f_superFA:s", 
	   "f_coil:s", 
	   "f_hamap:s", 
	   "f_pirsf:s", 
	  "splitXmlByID:s", # Used to prepare input of blast2go pipeline. 
	   "fitB2G!", 
	  "showV4convert!", 
	  "jnTSV_iprID!", # Used to join ipr-IDs and GO-IDs in in.iprV5.tsv file for futher usage. 
	   "fmt_raw!",    # In .raw format in iprV4, the columns are different. 
	   "srt_byChar!", # Sort annotations by ther character order. 
	   "rm_GOdesc!",  # Remove GO-description from iprV4 annotation if given. 
	# task="go"
	  "go_obo:s", # gene_ontology_edit_20150309.obo
	   "goIDColN:i", 
	   "addInfor:s", 
	   "go2ahrd:s", 
	"out:s", 
); 

sub usage {
	print <<HH; 
################################################################################
# perl $0 input_files
#
# -help 
#
# -task       [task_string] 
#              convert : 
#              -dbXml2GOlist : [Boolean] Convert interpro.xml to GO list 'GO:ID\\tCategory\\tDescription'
#               -xpath2go    : [/interprodb/interpro/class_list/classification] Normally no need to change. 
#
#              -xml5_to_raw4 : [Boolean] Convert interproV5.xml file to interproV4.raw file. 
#
#              -join_xml5    : [Boolean]
#
#              -showV4convert: [Boolean] Show the command to convert ipsV4_result from .raw to others. 
#               -dirV4new    : [/home/Sunhh/iprscan/iprscanV4_wiV5/] This is the new dir of iprV4 for converting files. 
#
#              -dbV5_to_dbV4 : [Boolean]
#               -dirV5       : [/share1/src/iprscan/iprscanV5/interproscan-5.11-51.0/] 
#               -dirV4old    : [/share1/src/iprscan/iprscanV4/iprscan/] This is the raw dir of iprV4 installation. 
#               -dirV4new    : [/home/Sunhh/iprscan/iprscanV4_wiV5/] This is the new dir of iprV4 for converting files. 
#               -hmmconvert  : [/share1/src/HMMER/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmconvert] Try to use the latest one because iprscan will update hmm files. 
#               -doIndex     : [Boolean] Run \$dirV4new/bin/index_data.pl to index the new database. 
#               -f_matXML    : [/share1/src/iprscan/iprscanV5/database_ipr_51.0/match_complete.xml.gz] 
#                               This should be downloaded from iprscan database. Could be with .gz
#               -f_matDTD    : [/share1/src/iprscan/iprscanV5/database_ipr_51.0/match_complete.dtd]
#               -f_iprXML    : [/share1/src/iprscan/iprscanV5/database_ipr_51.0/interpro.xml.gz] 
#                               This should be downloaded from iprscan database. Could be with .gz
#               -f_iprDTD    : [/share1/src/iprscan/iprscanV5/database_ipr_51.0/interpro.dtd]
#               -f_superF    : [\$dirV5/data/superfamily/1.75/hmmlib_1.75] 
#               -f_smrtTh    : [\$dirV5/data/smart/6.2/smart.thresholds]
#               -f_smrtDe    : [\$dirV5/data/smart/6.2/smart.desc]
#               -f_smrtHm    : [\$dirV5/data/smart/6.2/smart.HMMs]
#               -f_sfHMM     : [\$dirV5/data/pirsf/3.01/sf_hmm_all]
#               -f_prints    : [\$dirV5/data/prints/42.0/prints.pval]
#               -f_gene3d    : [\$dirV5/data/gene3d/3.5.0/gene3d_classified.hmm]
#               -f_tigrfam   : [\$dirV5/data/tigrfam/15.0/TIGRFAMs_15.0_HMM.LIB]
#               -f_pfamAH    : [\$dirV5/data/pfam/27.0/Pfam-A.hmm]
#               -f_pfamAS    : [\$dirV5/data/pfam/27.0/Pfam-A.seed]
#               -f_pfamC     : [\$dirV5/data/pfam/27.0/Pfam-C]
#               -f_prodom    : [\$dirV5/data/prodom/2006.1/prodom.ipr]
#               -f_panther   : [\$dirV5/data/panther/9.0/model/]
#               -f_prosEv    : [\$dirV5/data/prosite/20.105/evaluator.dat]
#               -f_prosite   : [\$dirV5/data/prosite/20.105/prosite.dat]
#               -f_fprint    : [\$dirV5/data/prints/42.0/FingerPRINTShierarchy.db]
#               -f_superFT   : [\$dirV5/data/superfamily/1.75/model.tab]
#               -f_superFA   : [\$dirV5/data/superfamily/1.75/model.tab]
#               -f_coil      : [\$dirV5/data/coils/2.2/new_coil.mat]
#               -f_hamap     : [\$dirV5/data/hamap/latest/hamap.prf]
#               -f_pirsf     : [\$dirV5/data/pirsf/3.01/pirsf.dat]
#
#              -splitXmlByID : [OutDir] Split ips_result.xml into small peices, in the way one ID.xml with one ID protein. 
#                               OutDir must not exist. 
#               -fitB2G      : [Boolean] Use this if you are splitting xml files for B2G. 
#                               This will add another node name 'EBIInterProScanResults'. 
# 
#              -jnTSV_iprID  : [Boolean] Merge ipr-IDs and GO-IDs in input.ipr.tsv file. 
#                              Ref : https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats 
#               -fmt_raw     : [Boolean] The input table is in iprV4.raw format instead of iprV5.tsv format. 
#               -srt_byChar  : [Boolean] Sort annotations by character order. 
#               -rm_GOdesc   : [Boolean] Remove GO description infor if given. (should work only for V4.raw)
#
#
#              go : 
#               -go_obo      : [gene_ontology_edit_20150309.obo] The GO .obo file being used. 
#               -goIDColN    : [0] The column number of GO ID. 
#               
#              -addInfor    : ['name'] Feature names being attached in lines. 
#
#              -go2ahrd     : [out_ahrd_pref] Create out_ahrd_pref.GOname and out_ahrd_pref.annot files for AHRD Terms: 
#                                                gene_ontology_result and blast2go 
# -out        [out_file_name] 
################################################################################
HH
	exit 1; 
}

# Basic settings. 
$opts{'help'} and &usage(); 
(keys %opts) == 0 and &usage(); 
our @InFp = (); 
if ( !@ARGV ) {
	@InFp = (\*STDIN); 
} else {
	push( @InFp, &openFH($_, '<') ) for (@ARGV); 
}

my $outFh = \*STDOUT; 
defined $opts{'out'} and $outFh = &openFH($opts{'out'}, '>'); 

my (%go_obo); 
if ( defined $opts{'go_obo'} ) {
	my $tmpFh = &openFH($opts{'go_obo'}, '<'); 
	my $tmp_class = 'Basic'; 
	my $tmp_ID = 'Basic'; 
	while (<$tmpFh>) {
		chomp; m/^\s*(#|$)/ and next; 
		if (m/^\[(\S+)\]$/) {
			$tmp_class = $1; 
			next; 
		}
		m/^([^\s:]+):\s+(\S.*)$/ or &stopErr("[Err] Bad GO.obo format : $_\n"); 
		my ($tag, $value)  = ($1, $2); 
		if ($tag =~ m/^id$/i) {
			$tmp_ID = $value; 
		}
		if ( $tag =~ m/^alt_id$/i ) {
			push(@{ $go_obo{'alt_id'}{$value} }, $tmp_ID); 
		}
		push(@{$go_obo{$tmp_class}{$tmp_ID}{$tag}}, $value); 
	}
	close ($tmpFh); 
	$opts{'goIDColN'} //= 0; 
}


$opts{'task'} //= ''; 

my %iprV4DbName; 
$iprV4DbName{'HAMAP'} = 'HAMAP';                        # (V4)HAMAP        (V5)Hamap                 : *High-quality Automated and Manual Annotation of Microbial Proteomes
$iprV4DbName{'PRODOM'} = 'blastprodom';                 # (V4)ProDom team  (V5)ProDom                : ProDom is a comprehensive set of protein domain families automatically generated from the UniProt Knowledge Database.
$iprV4DbName{'PROSITE_PROFILES'} = 'ProfileScan';       # (V4)PROSITE      (V5)ProSiteProfiles       : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them
$iprV4DbName{'PROSITE_PATTERNS'} = 'PatternScan';       # (V4)PROSITE      (V5)ProSitePatterns       : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them
$iprV4DbName{'PRINTS'} = 'FPrintScan';                  # (V4)PRINTS       (V5)PRINTS                : A fingerprint is a group of conserved motifs used to characterise a protein family
$iprV4DbName{'SUPERFAMILY'} = 'superfamily';            # (V4)SUPERFAMILY  (V5)SUPERFAMILY           : SUPERFAMILY is a database of structural and functional annotation for all proteins and genomes.
$iprV4DbName{'PANTHER'} = 'HMMPanther';                 # (V4)PANTHER      (V5)PANTHER               : The PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification System is a unique resource that classifies genes by their functions, using published scientific experimental evidence and evolutionary relationships to predict function even in the absence of direct experimental evidence.
$iprV4DbName{'GENE3D'} = 'Gene3D';                      # (V4)GENE3D       (V5)Gene3D                : Structural assignment for whole genes and genomes using the CATH domain structure database
$iprV4DbName{'PIRSF'} = 'HMMPIR';                       # (V4)PIR          (V5)PIRSF                 : The PIRSF concept is being used as a guiding principle to provide comprehensive and non-overlapping clustering of UniProtKB sequences into a hierarchical order to reflect their evolutionary relationships.
$iprV4DbName{'PFAM'} = 'HMMPfam';                       # (V4)Pfam         (V5)Pfam                  : A large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs)
$iprV4DbName{'TIGRFAM'} = 'HMMTigr';                    # (V4)TIGRFAMs     (V5)TIGRFAM               : TIGRFAMs are protein families based on Hidden Markov Models or HMMs
$iprV4DbName{'COILS'} = 'Coils';                        # (V4)Ncoils       (V5)Coils                 : Prediction of Coiled Coil Regions in Proteins
# Third party databases installed separately. 
# No Phobius databse result found.                      #                  (V5)Phobius               : A combined transmembrane topology and signal peptide predictor
$iprV4DbName{'PHOBIUS'} = 'PHOBIUS'; 
$iprV4DbName{'SMART'} = 'HMMSmart';                     # (V4)SMART (V5)SMART                 : SMART allows the identification and analysis of domain architectures based on Hidden Markov Models or HMMs
$iprV4DbName{'TMHMM'} = 'TMHMM';                        # (V4)TMHMM_v2.0 (V5)TMHMM                 : Prediction of transmembrane helices in proteins
$iprV4DbName{'SIGNALP_EUK'} = 'SignalP';                # (V4)SignalP_v3.0 (V5)SignalP_EUK           : SignalP (organism type eukaryotes) predicts the presence and location of signal peptide cleavage sites in amino acid sequences for eukaryotes.
$iprV4DbName{'SIGNALP_GRAM_POSITIVE'} = 'SignalP';      # (V4)SignalP_v3.0 (V5)SignalP_GRAM_POSITIVE : SignalP (organism type gram-positive prokaryotes) predicts the presence and location of signal peptide cleavage sites in amino acid sequences for gram-positive prokaryotes.
$iprV4DbName{'SIGNALP_GRAM_NEGATIVE'} = 'SignalP';      # (V4)SignalP_v3.0 (V5)SignalP_GRAM_NEGATIVE : SignalP (organism type gram-negative prokaryotes) predicts the presence and location of signal peptide cleavage sites in amino acid sequences for gram-negative prokaryotes.
# Not known database in V4 but absent in V5. 
$iprV4DbName{'Non_Seg'} = 'Seg';                        # (V4)Seg          (V5)?                     : Seg replaces low complexity regions in protein sequences with X characters. If a resulting protein sequence is used as a query for a BLAST search, the regions with X characters are ignored.
# There are also 'pfscan|scanregexp|seg' in interproscanV4 database types pointing to 'PROFILE|PROSITE|SEG', which can not be found in V5 data. But it doesn't matter because there has been 'PROSITE_PROFILES|PROSITE_PATTERNS' in V5. 
my %goType; 
$goType{'MOLECULAR_FUNCTION'} = 'Molecular Function'; 
$goType{'BIOLOGICAL_PROCESS'} = 'Biological Porcess'; 
$goType{'CELLULAR_COMPONENT'} = 'Cellular Component'; 

# Invoke sub-functions
&dbXml2GOlist() if ( $opts{'task'} eq 'convert' and $opts{'dbXml2GOlist'} ); 
&xml5_to_raw4() if ( $opts{'task'} eq 'convert' and $opts{'xml5_to_raw4'} ); 
&join_xml5() if ( $opts{'task'} eq 'convert' and $opts{'join_xml5'} ) ; 
&dbV5_to_dbV4() if ( $opts{'task'} eq 'convert' and $opts{'dbV5_to_dbV4'} ); 
&splitXmlByID() if ( $opts{'task'} eq 'convert' and defined $opts{'splitXmlByID'} ); 
&showV4convert() if ( $opts{'task'} eq 'convert' and $opts{'showV4convert'} ); 
&jnTSV_iprID() if ( $opts{'task'} eq 'convert' and defined $opts{'jnTSV_iprID'} ); 
&addInfor() if ( $opts{'task'} eq 'go' and defined $opts{'addInfor'} ); 
&go2ahrd() if ( $opts{'task'} eq 'go' and defined $opts{'go2ahrd'} ); 

# Sub-functions
sub join_xml5 {
	my $header_out = 0; 
	my $tail_txt = ''; 
	for my $fh ( @InFp ) {
		while (<$fh>) {
			if ( m!^\<\?xml version=! ) {
				my $to = $_; 
				my $tl = <$fh>; 
				$to .= $tl; 
				$tl =~ m!^\<protein\-matches! or &stopErr("[Err] Unknown format: $to"); 
				$header_out == 0 and print {$outFh} $to; 
				$header_out = 1; 
				next; 
			}
			if ( m!^\s*\</protein\-matches!o ) {
				$tail_txt = $_; 
				$tail_txt =~ m!^\</protein\-matches\>! or &stopErr("[Err] Unknown format: $tail_txt"); 
				next; 
			}
			print {$outFh} $_; 
		}
	}
	print {$outFh} $tail_txt; 
	close($outFh); 
	return; 
}# sub join_xml5 
sub jnTSV_iprID {
# Format of TSV : 
# https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats
# The TSV format presents the match data in columns as follows:
#
# C  1: Protein Accession (e.g. P51587)
# C  2: Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
# C  3: Sequence Length (e.g. 3418)
# C  4: Analysis (e.g. Pfam / PRINTS / Gene3D)
# C  5: Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
# C  6: Signature Description (e.g. BRCA2 repeat profile)
# C  7: Start location
# C  8: Stop location
# C  9: Score - is the e-value of the match reported by member database method (e.g. 3.1E-52)
# C 10: Status - is the status of the match (T: true)
# C 11: Date - is the date of the run
# C 12: (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprlookup option is switched on)
# C 13: (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprlookup option is switched on)
# C 14: (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
# C 15: (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
#
# Format of iprV4.raw : 
# ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/4/README.html
#
# C  1: NF00181542	is the id of the input sequence.
# C  2: 27A9BBAC0587AB84	is the crc64 (checksum) of the protein sequence (supposed to be unique).
# C  3: 272	is the length of the sequence (in AA).
# C  4: HMMPIR	is the anaysis method launched.
# C  5: PIRSF001424	is the database members entry for this match.
# C  6: Prephenate dehydratase	is the database member description for the entry.
# C  7: 1	is the start of the domain match.
# C  8: 270	is the end of the domain match.
# C  9: 6.5e-141	is the evalue of the match (reported by member database method).
# C 10: T	is the status of the match (T: true, ?: unknown).
# C 11: 06-Aug-2005	is the date of the run.
# C 12: IPR008237	is the corresponding InterPro entry (if iprlookup requested by the user).
# C 13: Prephenate dehydratase with ACT region	is the description of the InterPro entry.
# C 14: Molecular Function:prephenate dehydratase activity (GO:0004664)	is the GO (gene ontology) description for the InterPro entry.
#
#
	my %gene_infor; 
	my %cnt; 
	for my $fh ( @InFp ) {
		while (<$fh>) {
			$cnt{'line'}++; 
			m/^\s*(#|$)/ and next; 
			chomp; 
			my @ta = split(/\t/, $_); 
			$cnt{'gene_order'}{$ta[0]} //= $cnt{'line'}; 

			# For interpro-ID : Col = 11(ID) + 12(Description) 
			if (defined $ta[11] and $ta[11] ne '') {
				$ta[11] =~ m!^IPR\d+$! or &stopErr("[Err] Unknown iprID [$ta[11]] in line: $_\n"); 
				$gene_infor{$ta[0]}{ipr}{$ta[11]} //= [ $cnt{'line'}, "$ta[11] ($ta[12])" ]; 
			}
			unless ( defined $opts{'fmt_raw'} ) {
				# For interpro-ID : Col = 11(ID) + 12(Description) 
				# For GO-ID : Col = 13
				if (defined $ta[13] and $ta[13] ne '' and $ta[13] !~ m!^\s+$!) {
					for my $tb (split(/\s*\|\s*/, $ta[13])) {
						$tb =~ m!^GO:\d+$! or &stopErr("[Err] Unknown GO ID [$tb] in line: $_\n"); 
						$gene_infor{$ta[0]}{'go'}{$tb} //= [ $cnt{'line'}, "$tb ()" ]; 
					}
				}
				# For Pathways annotations : Col = 14
				if (defined $ta[14] and $ta[14] ne '' and $ta[14] !~ m!^\s+$!) {
					for my $tb (split(/\s*\|\s*/, $ta[14])) {
						$tb =~ m!^([^\s:]+):\s+\S+$! or &stopErr("[Err] Unknown PWY infor [$tb] in line: $_\n"); 
						$gene_infor{$ta[0]}{'pwy'}{$tb} //= [ $cnt{'line'}, $tb ]; 
					}
				}
			} else {
				# For interpro-ID : Col = 11(ID) + 12(Description) ; Same to V5; 
				# For GO-ID : Col = 13 
				if (defined $ta[13] and $ta[13] ne '' and $ta[13] !~ m!^\s+$!) {
					while ( $ta[13] =~ s!^(\S.*?)\s+\((GO:\d+)\),\s+!! ) {
						my ($desc, $id) = ($1, $2); 
						$opts{'rm_GOdesc'} and $desc = ''; 
						$gene_infor{$ta[0]}{'go'}{$id} //= [ $cnt{'line'}, "$id ($desc)" ]; 
					}
					if ($ta[13] ne '' and $ta[13] !~ m!^\s+$!) {
						$ta[13] =~ m!^(\S.*?)\s+\((GO:\d+)\)\s*$! or &stopErr("[Err] Bad format [$ta[13]] in line: $_\n"); 
						my ($desc, $id) = ($1, $2); 
						$opts{'rm_GOdesc'} and $desc = ''; 
						$gene_infor{$ta[0]}{'go'}{$id} //= [ $cnt{'line'}, "$id ($desc)" ]; 
					}
				}
				# For Pathways annotations : Col = NA 
			}
		}
	}

	print {$outFh} join("\t", qw/GeneID IPRs GOs PWYs/)."\n"; 
	for my $geneID ( sort { $cnt{'gene_order'}{$a} <=> $cnt{'gene_order'}{$b} } keys %{$cnt{'gene_order'}} ) {
		for my $tk (qw/ipr go pwy/) {
			$gene_infor{$geneID}{$tk} //= { '-fake' => [-1, ''] }; 
		}
		if (defined $opts{'srt_byChar'}) {
			print {$outFh} join("\t", 
			 $geneID, 
			 join( ', ', map { $gene_infor{$geneID}{'ipr'}{$_}[1] } sort { $a cmp $b } keys %{$gene_infor{$geneID}{'ipr'}} ), 
			 join( ', ', map { $gene_infor{$geneID}{'go' }{$_}[1] } sort { $a cmp $b } keys %{$gene_infor{$geneID}{'go' }} ), 
			 join( ';;', map { $gene_infor{$geneID}{'pwy'}{$_}[1] } sort { $a cmp $b } keys %{$gene_infor{$geneID}{'pwy'}} ) 
			)."\n"; 
		} else {
			print {$outFh} join("\t", 
			 $geneID, 
			 join( ', ', map { $gene_infor{$geneID}{'ipr'}{$_}[1] } sort { $gene_infor{$geneID}{'ipr'}{$a}[0] <=> $gene_infor{$geneID}{'ipr'}{$b}[0] } keys %{$gene_infor{$geneID}{'ipr'}} ), 
			 join( ', ', map { $gene_infor{$geneID}{'go' }{$_}[1] } sort { $gene_infor{$geneID}{'go' }{$a}[0] <=> $gene_infor{$geneID}{'go' }{$b}[0] } keys %{$gene_infor{$geneID}{'go' }} ), 
			 join( ';;', map { $gene_infor{$geneID}{'pwy'}{$_}[1] } sort { $gene_infor{$geneID}{'pwy'}{$a}[0] <=> $gene_infor{$geneID}{'pwy'}{$b}[0] } keys %{$gene_infor{$geneID}{'pwy'}} ) 
			)."\n"; 
		}
	}

	return(); 
}# jnTSV_iprID() 


sub go2ahrd {
	my $oname_Fh = &openFH("$opts{'go2ahrd'}.GOname", '>'); 
	my $oanno_Fh = &openFH("$opts{'go2ahrd'}.annot", '>'); 
	my %annotGO; 
	for my $fh ( @InFp ) {
		while (<$fh>) {
			chomp; m/^\s*(#|$)/ and next; 
			my @ta = split(/\t/, $_); 
			$ta[1] =~ m/^GO:/ or next; 
			if ( defined $ta[2] ) {
				if ( defined $annotGO{$ta[1]} ) {
					&tsmsg("[Wrn] Skip repeated annotation description of [$ta[1]]\n"); 
				} else {
					print {$oanno_Fh} join("\t", @ta[0,1,2])."\n"; 
				}
			}
			my $tr1 = &infor_by_goID($ta[1], \%go_obo); 
			if ( defined $tr1 ) {
				print {$oname_Fh} join("\t", $ta[0], 0.5, $ta[1], join(";; ", @{ $tr1->{'name'} }) )."\n"; 
			} else {
				print {$oname_Fh} join("\t", $ta[0], 0.5, $ta[1], 'NA')."\n"; 
			}
		}
		close($fh); 
	}
}# sub go2ahrd () 
sub addInfor {
	my @need_id = map { s!\s!!g; $_; } split(/,/, $opts{'addInfor'}); 
	for my $fh (@InFp) {
		while (<$fh>) {
			chomp; 
			my @ta = split(/\t/, $_); 
			my @added_cols; 
			my $go_id = $ta[ $opts{'goIDColN'} ]; 
			if ( $. == 1 and $go_id =~ m/^go(id)?$/i ) {
				print {$outFh} join("\t", $_, @need_id)."\n"; 
				next; 
			}
			my $tr1 = &infor_by_goID($go_id, \%go_obo); 
			if ( defined $tr1 ) {
				for my $tid (@need_id) {
					if ( defined $tr1->{$tid} ) {
						push(@added_cols, $tr1->{$tid}); 
					} else {
						push(@added_cols, ['']); 
					}
				}
			} else {
				for my $tid (@need_id) {
					push(@added_cols, ['']); 
				}
			}
			print {$outFh} "$_"; 
			for my $tr (@added_cols) {
				print {$outFh} "\t" . join(";; ", @$tr); 
			}
			print {$outFh} "\n"; 
		}
	}
}#sub addInfor () 
sub showV4convert {
	$opts{'dirV4new'} //= '/home/Sunhh/iprscan/iprscanV4_wiV5/'; 
	-d $opts{'dirV4new'} or &stopErr("[Err] Not valid interproscan V4.8 directory with -dirV4new\n"); 
	$opts{'dirV4new'} = &fileSunhh::_abs_path( $opts{'dirV4new'} ); 
	$opts{'dirV4new'} =~ s!/+$!!; 
	my $pl_cnvt = "$opts{'dirV4new'}/bin/"; 
}#sub showV4convert () 
sub splitXmlByID {
	$opts{'splitXmlByID'} //= 'OutDir'; 
	$opts{'splitXmlByID'} eq '' and $opts{'splitXmlByID'} = 'OutDir'; 
	-e $opts{'splitXmlByID'} and &stopErr("[Err] Already exists : |$opts{'splitXmlByID'}|\n"); 
	mkdir($opts{'splitXmlByID'}, 0755); 
	my $oDir = fileSunhh::_abs_path($opts{'splitXmlByID'}); 
	my %id2fn; 
	my $cnt = 0; 
	for my $fh (@InFp) {
		my ($header, $footer, $body); 
		my %status; 
		$status{'is_header'} = 1; 
		$status{'is_footer'} = 0; 
		$status{'is_body'}   = 0; 
		$status{'ID'} = []; 
		$opts{'fitB2G'} and $header = '<EBIInterProScanResults>' . "\n"; 
		while (<$fh>) {
			if ( $status{'is_header'} == 1 and $_ =~ m!^\s*\<interpro_matches\>\s*$|^\s*$|^\s*\<\?xml(\s+|\>)|^\s*\<protein\-matches(\s+|\>)! ) {
				$header .= $_; 
			} elsif ( ($status{'is_header'} == 1 or $status{'is_footer'} == 1 or $status{'is_body'} == 0) and $_ =~ m!^\s*\<protein(\>|\s+)! ) {
				$status{'is_body'} = 1; 
				$status{'is_footer'} = 0; 
				$status{'is_header'} = 0; 
				$body = $_; 
				if ( $_ =~ m!\sid=\"([^"\s]+)\"!i ) {
					push( @{$status{'ID'}} , $1 ); 
				}
			} elsif ( $status{'is_body'} == 1 and $_ =~ m!^\s*\<\/protein\>! ) {
				$body .= $_; 
				@{$status{'ID'}} == 0 and &tsmsg("[Wrn] No protein ID captured in body:\n$body\n"); 
				grep { defined $id2fn{ $_ } and &stopErr("[Err] repeated ID $_ to output.\n") } @{$status{'ID'}} ; 
				for my $id ( @{$status{'ID'}} ) {
					$id eq '' and next; 
					my $ofn = "$oDir/${id}.xml"; 
					$id2fn{$id} = $ofn; 
					open O,'>',"$ofn" or &stopErr("[Err] Failed to write $ofn\n"); 
					print O "${header}${body}"; 
					close O; 
				}
				$status{'is_body'} = 0; 
				$body = ''; 
				$status{'ID'} = []; 
			} elsif ( $status{'is_body'} == 0 and $_ =~ m/^\s*$|^\s*\<\/protein\-matches(\s+|\>)/ ) {
				$footer = $_; 
				$status{'is_footer'} = 1; 
			} elsif ( $status{'is_body'} == 1 ) {
				if ( m!^\s*\<xref\s+id=\"([^"\s]+)\"!i ) {
					push(@{$status{'ID'}}, $1); 
				}
				$body .= $_; 
			} elsif ( $status{'is_footer'} == 1 ) {
				$footer .= $_; 
			} else {
				&stopErr("[Err] Unknown line: $_\n"); 
			}
		}
		close($fh); 
		$opts{'fitB2G'} and $footer .= "</EBIInterProScanResults>\n"; 
		for my $id (sort keys %id2fn) {
			open OO,'>>',"$id2fn{$id}" or &stopErr("[Err] Failed to write $id2fn{$id}\n"); 
			print OO $footer; 
			close OO; 
		}
	}
}# End sub splitXmlByID() 

sub dbV5_to_dbV4 {
	$opts{'dirV5'} //= '/share1/src/iprscan/iprscanV5/interproscan-5.11-51.0/'; 
	$opts{'dirV4old'} //= '/share1/src/iprscan/iprscanV4/iprscan/'; 
	$opts{'dirV4new'} //= '/home/Sunhh/iprscan/iprscanV4_wiV5/'; 
	-e $opts{'dirV4new'} and &stopErr("[Err] Already exists $opts{'dirV4new'}\n"); 
	mkdir($opts{'dirV4new'}, 0755); 
	$opts{$_} = fileSunhh::_abs_path($opts{$_}) foreach ( qw/dirV5 dirV4old dirV4new/ ); 
	$opts{'dirV5'} =~ s!/+$!!; 
	$opts{'dirV4old'} =~ s!/+$!!; 
	$opts{'dirV4new'} =~ s!/+$!!; 
	
	$opts{'hmmconvert'} //= '/share1/src/HMMER/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmconvert'; 
	$opts{'f_matXML'} //= '/share1/src/iprscan/iprscanV5/database_ipr_51.0/match_complete.xml.gz'; 
	$opts{'f_matDTD'} //= '/share1/src/iprscan/iprscanV5/database_ipr_51.0/match_complete.dtd'; 
	$opts{'f_iprXML'} //= '/share1/src/iprscan/iprscanV5/database_ipr_51.0/interpro.xml.gz'; 
	$opts{'f_iprDTD'} //= '/share1/src/iprscan/iprscanV5/database_ipr_51.0/interpro.dtd'; 
	$opts{'f_superF'} //= $opts{'dirV5'} . '/data/superfamily/1.75/hmmlib_1.75'; 
	$opts{'f_smrtTh'} //= $opts{'dirV5'} . '/data/smart/6.2/smart.thresholds'; 
	$opts{'f_smrtDe'} //= $opts{'dirV5'} . '/data/smart/6.2/smart.desc'; 
	$opts{'f_smrtHm'} //= $opts{'dirV5'} . '/data/smart/6.2/smart.HMMs'; 
	$opts{'f_sfHMM'}  //= $opts{'dirV5'} . '/data/pirsf/3.01/sf_hmm_all'; 
	$opts{'f_prints'} //= $opts{'dirV5'} . '/data/prints/42.0/prints.pval'; 
	$opts{'f_gene3d'} //= $opts{'dirV5'} . '/data/gene3d/3.5.0/gene3d_classified.hmm'; 
	$opts{'f_tigrfam'} //= $opts{'dirV5'} . '/data/tigrfam/15.0/TIGRFAMs_15.0_HMM.LIB'; 
	$opts{'f_pfamAH'} //= $opts{'dirV5'} . '/data/pfam/27.0/Pfam-A.hmm'; 
	$opts{'f_pfamAS'} //= $opts{'dirV5'} . '/data/pfam/27.0/Pfam-A.seed'; 
	$opts{'f_pfamC'}  //= $opts{'dirV5'} . '/data/pfam/27.0/Pfam-C'; 
	$opts{'f_prodom'} //= $opts{'dirV5'} . '/data/prodom/2006.1/prodom.ipr'; 
	$opts{'f_panther'} //= $opts{'dirV5'} . '/data/panther/9.0/model/'; 
	$opts{'f_prosEv'}  //= $opts{'dirV5'} . '/data/prosite/20.105/evaluator.dat'; 
	$opts{'f_prosite'} //= $opts{'dirV5'} . '/data/prosite/20.105/prosite.dat'; 
	$opts{'f_fprint'}  //= $opts{'dirV5'} . '/data/prints/42.0/FingerPRINTShierarchy.db'; 
	$opts{'f_superFT'} //= $opts{'dirV5'} . '/data/superfamily/1.75/model.tab'; 
	$opts{'f_superFA'} //= $opts{'dirV5'} . '/data/superfamily/1.75/model.tab'; 
	$opts{'f_coil'}    //= $opts{'dirV5'} . '/data/coils/2.2/new_coil.mat'; 
	$opts{'f_hamap'}   //= $opts{'dirV5'} . '/data/hamap/latest/hamap.prf'; 
	$opts{'f_pirsf'}   //= $opts{'dirV5'} . 'data/pirsf/3.01/pirsf.dat'; 

	# Construct dirV4new . 
	mkdir("$opts{'dirV4new'}/data/", 0755); 
	mkdir("$opts{'dirV4new'}/tmp/", 0777); 
	mkdir("$opts{'dirV4new'}/bin/", 0755); 
	mkdir("$opts{'dirV4new'}/conf/", 0755); 
	for my $name ( qw/Config.pl images lib/ ) {
		symlink("$opts{'dirV4old'}/$name", "$opts{'dirV4new'}/$name") or &stopErr("[Err] Failed to make link for $name\n"); 
	}
	my %binFiles = map { $_=>1 } qw/converter.pl gene3d.pl index_data.pl iprscan iprscan_wrapper.pl meter.pl out2raw.pl pantherScore.pl ProDomBlast3i.pl ResubmitJobs.pl wget.pl/; 
	opendir DD,"$opts{'dirV4old'}/bin/" or &stopErr("[Err] failed to read $opts{'dirV4old'}/bin/\n"); 
	for my $fn ( grep { $_ !~ m!^\.! } readdir(DD) ) {
		defined $binFiles{$fn} or do { symlink("$opts{'dirV4old'}/bin/$fn", "$opts{'dirV4new'}/bin/$fn") or &stopErr("[Err] Failed to link file $fn\n"); next; }; 
		open F,'<',"$opts{'dirV4old'}/bin/$fn" or &stopErr("[Err] Failed to read file $fn\n"); 
		open O,'>',"$opts{'dirV4new'}/bin/$fn" or &stopErr("[Err] Failed to write file $fn\n"); 
		while (<F>) {
			if ( $_ =~ m!^(\s*\$ENV\{\s*['"]?IPRSCAN_HOME['"]?\s*\}\s*=\s*['"])[^\s'"]*(['"]\s*;.*)$! ) {
				$_ = "$1$opts{'dirV4new'}$2"; 
				chomp($_); $_ .= "\n"; 
			}
			print O $_; 
		}
		close O; 
		close F; 
	}
	closedir DD; 
	my %confFiles = map { $_=>1 } qw/sixpack.sh seqret.sh/; 
	opendir DD,"$opts{'dirV4old'}/conf/" or &stopErr("[Err] failed to read conf/\n"); 
	for my $fn ( grep { $_ !~ m!^\.! } readdir(DD) ) {
		defined $confFiles{$fn} or do { symlink("$opts{'dirV4old'}/conf/$fn", "$opts{'dirV4new'}/conf/$fn") or &stopErr("[Err] Failed to link file $fn\n"); next; }; 
		open F,'<',"$opts{'dirV4old'}/conf/$fn" or &stopErr("[Err] Failed to read file $fn\n"); 
		open O,'>',"$opts{'dirV4new'}/conf/$fn" or &stopErr("[Err] Failed to write file $fn\n"); 
		while (<F>) {
			if ( $_ =~ m!^(\s*IPRSCAN_HOME\s*=\s*).+?$! ) {
				$_ = "$1$opts{'dirV4new'}"; 
				chomp($_); $_ .= "\n"; 
			}
			print O $_; 
		}
		close O; 
		close F; 
	}
	closedir DD; 
	# Generate data/ in $dirV4new
	my %out; 
	$out{'f_matXML'}  //= 'match_complete.xml';
	$out{'f_matDTD'}  //= 'match_complete.dtd';
	$out{'f_iprXML'}  //= 'interpro.xml';
	$out{'f_iprDTD'}  //= 'interpro.dtd';
	$out{'f_superF'}  //= 'superfamily.hmm'; 
	$out{'f_smrtTh'}  //= 'smart.thresholds'; 
	$out{'f_smrtDe'}  //= 'smart.desc'; 
	$out{'f_smrtHm'}  //= 'smart.HMMs'; 
	$out{'f_sfHMM'}   //= 'sf_hmm'; 
	$out{'f_prints'}  //= 'prints.pval'; 
	$out{'f_gene3d'}  //= 'gene3d.lib'; 
	$out{'f_tigrfam'} //= 'TIGRFAMs_HMM.LIB'; 
	$out{'f_pfamAH'}  //= 'Pfam-A.hmm'; 
	$out{'f_pfamAS'}  //= 'Pfam-A.seed'; 
	$out{'f_pfamC'}   //= 'Pfam-C'; 
	$out{'f_prodom'}  //= 'prodom.ipr'; 
	$out{'f_panther'} //= 'Panther'; 
	$out{'f_prosEv'}  //= 'evaluator.dat';
	$out{'f_prosite'} //= 'prosite.dat';
	$out{'f_fprint'}  //= 'FingerPRINTShierarchy.db';
	$out{'f_superFT'} //= 'superfamily.tab';
	$out{'f_superFA'} //= 'superfamily.acc';
	$out{'f_coil'}    //= 'new_coil.mat';
	$out{'f_hamap'}   //= 'hamap.prf'; 
	$out{'f_pirsf'}   //= 'pirsf.dat'; 
	$out{$_} = "$opts{'dirV4new'}/data/$out{$_}" foreach ( keys %out ); 
	# Still we have following files not found: sf.seq sf.tb sf_hmm_sub model2sf_map.csv superfamily.acc
	# While superfamily.acc is the 1st column of superfamily.tab, and the other files should have no changes. 
	for my $fn (qw/sf.seq sf.tb sf_hmm_sub model2sf_map.csv/) {
		symlink("$opts{'dirV4old'}/data/$fn", "$opts{'dirV4new'}/data/$fn"); 
	}
	my %bin_hmm3 = map { $_=>1 } qw/f_pfamAH f_tigrfam f_gene3d/; 
	my %bin_hmm2 = map { $_=>1 } qw/f_sfHMM f_superF f_sfHMM/; 
	for my $fnk ( keys %out ) {
		defined $opts{$fnk} or &stopErr("[Err] No $fnk ($out{$fnk}) source defined.\n"); 
		my ( $fn, $path_aref ) = &fileSunhh::_chkExist($opts{$fnk}); 
		( defined $fn and $fn ne '' ) or do { &tsmsg("[Err] No source file: $opts{$fnk}\n"); next; }; 
		if ( $fnk eq 'f_superFA' ) {
			&exeCmd_1cmd("cut -f 1 $opts{$fnk} > $out{$fnk}"); 
			next; 
		} elsif ( $fnk =~ m!^(f_matXML|f_matDTD|f_iprXML|f_iprDTD)$! ) {
			if ( $opts{$fnk} =~ m!\.gz(ip)?\s*$!i ) {
				&exeCmd_1cmd("gzip -cd $opts{$fnk} > $out{$fnk}"); 
			} else {
				symlink( $opts{$fnk}, $out{$fnk} ) or &stopErr("[Err] Failed to link $opts{$fnk} to $out{$fnk}\n"); 
			}
			next; 
		}
		if ( -B $opts{$fnk} ) {
			symlink( $opts{$fnk}, $out{$fnk} ) or &stopErr("[Err] Failed to link $opts{$fnk} to $out{$fnk}\n"); 
			next; 
		}
		if ( defined $bin_hmm3{$fnk} ) {
			open F,'<',"$opts{$fnk}" or &stopErr("[Err] Failed to read $fnk\n"); 
			my $tmp_hh = <F>; 
			close F; 
			$tmp_hh =~ m!^HMMER3\b!i or &stopErr("[Err] Head line of $opts{$fnk} format wrong\n"); 
			symlink( $opts{$fnk}, $out{$fnk} ) or &stopErr("[Err] Failed to link $opts{$fnk} to $out{$fnk}\n"); 
		} elsif ( defined $bin_hmm2{$fnk} ) {
			open F,'<',"$opts{$fnk}" or &stopErr("[Err] Failed to read $fnk\n"); 
			my $tmp_hh = <F>; 
			close F; 
			if ( $tmp_hh =~ m!^HMMER2\b!i ) {
				symlink( $opts{$fnk}, $out{$fnk} ) or &stopErr("[Err] Failed to link $opts{$fnk} to $out{$fnk}\n"); 
			} else {
				&exeCmd_1cmd("$opts{'hmmconvert'} -2 $opts{$fnk} > $out{$fnk}"); 
			}
		} else {
			symlink( $opts{$fnk}, $out{$fnk} ) or &stopErr("[Err] Failed to link $opts{$fnk} to $out{$fnk}\n"); 
		}
	}#End for my $fnk 
	# Run index_data.pl
	if ( $opts{'doIndex'} ) {
		&exeCmd_1cmd("perl $opts{'dirV4new'}/bin/index_data.pl -inx -bin -v -p $opts{'dirV4new'}/data/"); 
	}
}# sub dbV5_to_dbV4 () 



sub dbXml2GOlist {
	$opts{'xpath2go'} //= undef(); 
	my %goInfor = %{ &goFromXML_fh(@InFp, 'xpath2go'=>$opts{'xpath2go'}) }; 
	print {$outFh} join("\t", qw/GO_ID Category Description/)."\n"; 
	for my $id ( sort keys %goInfor ) {
		print {$outFh} join("\t", $id, $goInfor{$id}{'category'}, $goInfor{$id}{'description'})."\n"; 
	}
}# sub dbXml2GOlist() 

# It seems the ipsV4 use subject_ID to trace back database, so if the subject_ID doesn't exist in ipsV4 database, 
#  there will be problem. So a better way is to write iprV4_result.xml directly. 
sub xml5_to_raw4 {
	# my ($parm_href, @fileFH) = &_parmFromFH( @_ ); 
	# my %parm = %$parm_href; undef($parm_href); 
	my $parser = XML::LibXML->new(); 
	for my $fh (@InFp) {
		my $doc = $parser->parse_fh($fh); 
		for my $node ( $doc->findnodes('*') ) {
			for my $prot_node ( $node->findnodes('*') ) {
				$prot_node->nodeName() eq 'protein' or &stopErr("[Err] <protein> is not the 2nd tier of xml.\n"); 
				my %tag2val; 
				for my $child_node ( $prot_node->findnodes('*') ) {
					my $child_name = $child_node->nodeName(); 
					if ( $child_name eq 'sequence' ) {
						$tag2val{'md5'} = $child_node->getAttribute('md5'); 
						$tag2val{'seq'} = $child_node->textContent(); 
						$tag2val{'len'} = length( $tag2val{'seq'} ); 
					} elsif ( $child_name eq 'xref' ) {
						push(@{$tag2val{'id'}}, $child_node->getAttribute('id')); 
						push(@{$tag2val{'name'}}, $child_node->getAttribute('name')); 
					} elsif ( $child_name eq 'matches' ) {
						for my $match_node ( $child_node->findnodes('*') ) {
							my $mat2href = &_match2Hash( $match_node ); 
							my $rawAref  = &_matchHref2rawAref( $mat2href ); 
							$rawAref->[0] = $tag2val{'id'}; 
							$rawAref->[1] = $tag2val{'md5'}; 
							$rawAref->[2] = $tag2val{'len'}; 
							defined $iprV4DbName{ $rawAref->[3] } or do { &tsmsg( "[Wrn] Undefined DB name $rawAref->[3]\n" ); $iprV4DbName{ $rawAref->[3] } = $rawAref->[3];  }; 
							$rawAref->[3] = $iprV4DbName{ $rawAref->[3] } // $rawAref->[3]; 
							$rawAref->[10] = '15-4-2015'; 
							for my $tmpID ( @{$rawAref->[0]} ) {
								for (my $i=0; $i < @{ $rawAref->[6] }; $i++ ) {
									print {$outFh} join("\t", $tmpID, @{$rawAref}[1..5], 
									 $rawAref->[6][$i], 
									 $rawAref->[7][$i], 
									 $rawAref->[8][$i], 
									 @{$rawAref}[9..13]
									)."\n"; 
								}
							}#End for my $tmpID
						}
					}
				}#End for my $child_node 
			}
		}#End for my $node 
	}
}# sub xml5_to_raw4() 

sub xml5_to_xml4 {
	# Nothing happens yet. 
	# I want to invoke iprscanV5 scripts and iprscanV4 scripts to finish this work. 
	# Step1. I need to convert iprscanV5_result.xml to iprscanV5_result.tsv 
	# Step2. Use iprscanV4 to convert iprscanV5_result.tsv to iprscanV5_result.iprV4.xml
	# Between Step1 and Step2, I want to re-index iprscanV4 databases with iprscanV5 databases, 
	#  in which way I think it will help convert more annotations. 
}# sub xml5_to_xml4() 




=head1 _matchHref2rawAref( &_match2Hash($in_match_node) )

ipsV4_raw format: 
InterProScan makes results available in four formats {raw ebixml xml txt html}:

* raw format
       - is basic tab delimited format useful for uploading the data into a
         relational database or concatenation of different runs.
       - is all on one line.
       - Example here (with descriptions):
--------------------------------------------------------------------------------
 NF00181542      0A5FDCE74AB7C3AD        272     HMMPIR  PIRSF001424     Prephenate dehydratase  1       270     6.5e-141        T       06-Aug-2005         IPR008237       Prephenate dehydratase with ACT region  Molecular Function:prephenate dehydratase activity (GO:0004664), Biological Process:L-phenylalanine biosynthesis (GO:0009094)
 
       Where: NF00181542:             is the id of the input sequence.
              27A9BBAC0587AB84:       is the crc64 (checksum) of the protein sequence (supposed to be unique).
              272:                    is the length of the sequence (in AA).
              HMMPIR:                 is the anaysis method launched.
              PIRSF001424:            is the database members entry for this match.
              Prephenate dehydratase: is the database member description for the entry.
              1:                      is the start of the domain match.
              270:                    is the end of the domain match.
              6.5e-141:               is the evalue of the match (reported by member database method).
              T:                      is the status of the match (T: true, ?: unknown).
              06-Aug-2005:            is the date of the run.
              IPR008237:              is the corresponding InterPro entry (if iprlookup requested by the user).
              Prephenate dehydratase with ACT region:                           is the description of the InterPro entry.
              Molecular Function:prephenate dehydratase activity (GO:0004664):  is the GO (gene ontology) description for the InterPro entry.
--------------------------------------------------------------------------------
Return    : (\@back_arr)
 0-2   : ''; 
 3-5   : String // ''; 
 6-8   : [ [Start, End, Evalue], [], ... ]
 9-10  : 'T', ''; 
 11-12 : IPR_ID, IPR_description; 
 13    : merged GO annot: "GO_category:GO_description (GO:IDnum)", "Section1, Section2, ..."; 
 14    : pathway infor if any. "KEGG: 00230+2.7.7.6|KEGG: 00240+2.7.7.6|Reactome: REACT_1788". This is not included in interproV4 .raw result. 
=cut
sub _matchHref2rawAref {
	my $href = shift; 
	my %h = %$href; 
	my @back_arr; 
	$back_arr[0] = ''; # Query ID not here. 
	$back_arr[1] = ''; # md5/crc64 not here. 
	$back_arr[2] =  0; # Query sequence length not here. 
	$back_arr[3] = $h{'db'}{'library'} // &stopErr("[Err] No anaysis method found.\n"); 
	$back_arr[4] = $h{'db'}{'ac'} // &stopErr("[Err] No database members entry found.\n"); 
	$back_arr[5] = $h{'db'}{'desc'} // ''; 
	$back_arr[9] = 'T'; 
	$back_arr[10] = ''; 
	for my $th1 ( @{$h{'loc'}} ) {
		push(@{$back_arr[6]}, $th1->{'start'} // ''); 
		push(@{$back_arr[7]}, $th1->{'end'}   // ''); 
		push(@{$back_arr[8]}, $th1->{'evalue'} // $h{'evalue'} // ''); 
		# push(@{$back_arr[10]}, 'T'); # This is always "T"
		# push(@{$back_arr[11]}, ''}); # This is the date of input file. Value is localtime((stat $file)[9])->dmy ; 
	}
	if ( !( defined $back_arr[6] ) ) {
		$back_arr[6] = ['']; 
		$back_arr[7] = ['']; 
		$back_arr[8] = [ $h{'evalue'} // '']; 
	}
	if ( defined $h{'ipr'} ) {
		$back_arr[11] = $h{'ipr'}{'ac'}; 
		$back_arr[12] = $h{'ipr'}{'desc'}; 
		$back_arr[13] = ''; 
		if ( defined $h{'go'} ) {
			my @go_arr_1; 
			# my @go_arr_2; 
			for my $th2 ( @{$h{'go'}} ) {
				my $go_type = $goType{$th2->{'category'}} // $th2->{'category'}; 
				push(@go_arr_1, "$go_type:$th2->{'name'} ($th2->{'id'})"); 
				# push(@go_arr_2, {%$th2}); 
			}
			$back_arr[13] = join(', ', @go_arr_1); 
			# $back_arr[13][0] = join(', ', @go_arr_1); 
			# $back_arr[13][1] = \@go_arr_1; 
			# $back_arr[13][2] = \@go_arr_2; 
		} else {
			$back_arr[13] = ''; 
			# $back_arr[13][0] = ''; 
		}
	} else {
		$back_arr[11] = ''; 
		$back_arr[12] = ''; 
		$back_arr[13] = ''; 
	}
	if ( defined $h{'pathway'} ) {
		my @pwy_arr; 
		for my $th3 ( sort { $a->{'id'} cmp $b->{'id'} } @{$h{'pathway'}} ) {
			push(@pwy_arr, "$th3->{'db'}: $th3->{'id'}"); 
		}
		$back_arr[14] = join('|', @pwy_arr); 
	} else {
		$back_arr[14] = ''; 
	}

	return (\@back_arr); 
}# sub _matchHref2rawAref() 

# Supporting functions 
=head1 _match2Hash($match_node) # From <hmmer3-match|hmmer2-match|panther-match ...>, child of <matches>

Input      : XML::LibXML::Node object for child of <matches> from interproV5.xml
Return     : \%tag2val
 {'method'} = $match_node->nodeName(); # Align method. 
 {'evalue'} = the evalue of the match; # Not in 'COILS|PHOBIUS|SIGNALP_EUK|SIGNALP_GRAM_NEGATIVE|SIGNALP_GRAM_POSITIVE|TMHMM'
 {'score'}  = the score  of the match; # Not in 'COILS|PHOBIUS|SIGNALP_EUK|SIGNALP_GRAM_NEGATIVE|SIGNALP_GRAM_POSITIVE|TMHMM'
 Common {'db'}  => {'library|version|ac|name|(desc)?'}
 Common {'ipr'} => {'ac|desc|name'}, 
 Common  {'go'} => [ {'category|db|id|name'}=>Values ]
 Common {'pathway'} => [ {'db|id|name'} ]
 Common {'loc'} => [ {'start|end|(score|evalue)?'}, {}, {}]
 When {'db'}{'library'} : 
   == 'SMART' : {'db'} => {'ac|desc|name'} ; {'loc'} - {                  'hmm-start|hmm-end|hmm-length'}
   == 'GENE3D': {'db'} => {'ac'} ;           {'loc'} - {'env-end|env-start|hmm-start|hmm-end|hmm-length'}
   == 'PFAM'  : {'db'} => {'ac|desc|name'} ; {'loc'} - {'env-end|env-start|hmm-start|hmm-end|hmm-length'}
   == 'PANTHER':{'db'} => {'ac|     name'} ; {'loc'} - {''}
   == 'SUPERFAMILY': {'db'}=>{'ac|  name'} ; {'loc'} - {''}
   == 'PROSITE_PATTERNS': {'db'}=>{'ac|desc|name'}; {'loc'} - {'level'} with Child <alignment>
   == 'PROSITE_PROFILES': {'db'}; {'loc'} - {''} with Child <alignment>
   == 'COILS': {'db'} 
   == 'PHOBIUS': 
   == 'SIGNALP_EUK': {'loc'} with 'score'
   == 'SIGNALP_GRAM_NEGATIVE': with 'score'
   == 'SIGNALP_GRAM_POSITIVE': with 'score'
   == 'TMHMM'
   == 'TIGRFAM': {'loc'} - {'env-end|env-start|hmm-start|hmm-end|hmm-length'}
   == 'PRINTS':  {'loc'} - {'motifNumber|pvalue'} could be multiple children. 
   == 'HAMAP' :  {'loc'} with Child <alignment>
   == 'PIRSF' :  {'loc'} - {'env-end|env-start|hmm-start|hmm-end|hmm-length'}
   == 'PRODOM':  {'loc'}

=cut
sub _match2Hash {
	my $match_node = shift; 
	my %tag2val; 
	$tag2val{'method'} = $match_node->nodeName(); # <hmmer3-match| ...>
	$tag2val{'evalue'} = $match_node->getAttribute('evalue'); 
	$tag2val{'score'}  = $match_node->getAttribute('score'); 
	for my $child_node ( $match_node->findnodes('*') ) {
		# $child_node == <signature|locations ...>
		if ( $child_node->nodeName() eq 'signature' ) {
			for my $attr ( map { $_->localname() } $child_node->attributes('') ) {
				$tag2val{'db'}{$attr} = $child_node->getAttribute($attr); # ac/desc/name
			}
			for my $node2 ( $child_node->findnodes('*') ) {
				# $node2 == <entry|models|signature-library-release ...>
				if ( $node2->nodeName() eq 'entry' ) {
					for my $attr ( map { $_->localname() } $node2->attributes() ) {
						$tag2val{'ipr'}{$attr} = $node2->getAttribute($attr); # ac/desc/name/type
					}
					for my $node3 ( $node2->findnodes('*') ) {
						if ( $node3->nodeName() eq 'go-xref' ) {
							push(@{$tag2val{'go'}}, {}); 
							for my $attr ( map { $_->localname() } $node3->attributes() ) {
								$tag2val{'go'}[-1]{$attr} = $node3->getAttribute($attr); 
							}
						} elsif ( $node3->nodeName() eq 'pathway-xref' ) {
							push(@{$tag2val{'pathway'}}, {}); 
							for my $attr ( map { $_->localname() } $node3->attributes() ) {
								$tag2val{'pathway'}[-1]{$attr} = $node3->getAttribute($attr); 
							}
						} else {
							&stopErr("[Err] nodeName ", $node3->nodeName()," is not go-xref\n"); 
						}
					}
				} elsif ( $node2->nodeName() eq 'models' ) {
					; 
				} elsif ( $node2->nodeName() eq 'signature-library-release' ) {
					for my $attr ( map { $_->localname() } $node2->attributes() ) {
						$tag2val{'db'}{$attr} = $node2->getAttribute($attr); # library/version
					}
				} else {
					&stopErr("[Err] nodeName=",$node2->nodeName(), "\n"); 
				}
				# $node2 == </entry|models|signature-library-release>
			}
		} elsif ( $child_node->nodeName() eq 'locations' ) {
			for my $node2 ( $child_node->findnodes('*') ) {
				push(@{$tag2val{'loc'}}, {}); 
				for my $attr ( map { $_->localname() } $node2->attributes() ) {
					$tag2val{'loc'}[-1]{$attr} = $node2->getAttribute($attr); 
				}
			}
		} else {
			&stopErr("[Err] ", $child_node->nodeName(), " nodeName unknown.\n"); 
		}
		# $child_node == </signature|locations>
	}# End for my $child_nodue
	
	return (\%tag2val); 
}# sub _match2Hash () 



# Input    : ($inFh1, $inFh2, 'xpath2go'=>'/interprodb/interpro/class_list/classification')
# Return   : ( \%goInfor )
#  {$GO_ID}{'category'} => 'Molecular Function' 
#  {$GO_ID}{'description'} => 'DNA binding'
sub goFromXML_fh {
	my ($parm_href, @fileFH) = &_parmFromFH( @_ ); 
	my %parm = %$parm_href; undef($parm_href); 
	$parm{'xpath2go'} //= '/interprodb/interpro/class_list/classification'; 
	my $parser = XML::LibXML->new(); 
	my %has_go; 
	for my $fh ( @fileFH ) {
		my $doc = $parser->parse_fh($fh); 
		for my $node ( $doc->findnodes($parm{'xpath2go'}) ) {
			my $id = $node->getAttribute("id"); 
			defined $has_go{$id} and next; 
			for my $child ( $node->findnodes('*') ) {
				$has_go{$id}{ $child->nodeName() } = $child->textContent(); 
			}
		}
	}
	return (\%has_go); 
}# sub goFromXML_fh () 

sub _parmFromFH {
	my @fh; 
	my %parm; 
	for (my $i=0; $i<@_; $i++) {
		if ( ref($_[$i]) eq '' ) {
			$parm{$_[$i]} = $_[$i+1]; 
			$i ++; 
		} else {
			ref($_[$i]) eq 'GLOB' or &stopErr("[Err] [$_[$i]] is not a file handle\n"); 
			push(@fh, $_[$i]); 
		}
	}
	return(\%parm, @fh); 
}# sub _parmFromFH() 

# Input   : (\%hash_ref_1, \%hash_ref_2)
# Return  : (\%joined_hash)
sub joinHash {
	my %back; 
	for my $th ( @_ ) {
		for ( keys %$th ) {
			defined $back{$_} and next; 
			$back{$_} = $th->{$_}; 
		}
	}
	return (\%back); 
}# joinHash() 

=head1 infor_by_goID($go_ID, \%go_obo)

Return        : $go_obo{'Term'}{$go_ID} / undef() 
=cut
sub infor_by_goID {
	my ($goID, $go_obo_href) = @_; 
	$go_obo_href //= \%go_obo; 
	my $back_href = undef(); 
	$goID = $go_obo_href->{'alt_id'}{$goID}[0] // $goID; 
	if ( defined $go_obo_href->{'Term'}{$goID} ) {
		$back_href = $go_obo_href->{'Term'}{$goID}; 
	}
	return $back_href; 
}# sub infor_by_goID() 
