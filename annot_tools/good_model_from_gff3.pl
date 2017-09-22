#!/usr/bin/perl
# Assume you have filtered incomplete proteins and redundant proteins. 
# Now we need to filter genes too near to another gene/models. 
# I want to use this script as a shell. 
# 20150411 Use this script to get good gene models for training ab initio predictors. 
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"seedGff:s", # The gff3 file with good gene models satisfying complete and non-redundant rules. 
	"bgGff:s@", # background gff3 files, which is used to find models as islands. 
	"flankSize:i@", # flanking sizes which is used to find models as islands
	"flankName:s@", # Tag name of different flanking sizes. Will use flanking size if not given. 
	"pl_dg:s", # Path to deal_gff3.pl script. 
); 
sub usage {
	print <<HH;
###############################################################################
# perl $0 -seedGff prev_good.gff3 
#  -help 
#  This is a temporary script, and I will write a better one later. 
#  
#  -bgGff        [background.gff3] could be multiple times. 
#  -flankSize    [2000] could be multiple times. 
#  -flankName    [SameToFlankSize] could be multiple times. 
#  
#  -pl_dg        [/home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl]
#                  Path to deal_gff3.pl 
###############################################################################
HH
	exit 1; 
}

$opts{'help'} and &usage(); 
defined $opts{'seedGff'} or &usage(); 
my (@bgGffs, @fl_sizes, @fl_names); 
@bgGffs = @{$opts{'bgGff'}}; 
@fl_sizes = (defined $opts{'flankSize'}) ? @{$opts{'flankSize'}} : (2000) ; 
{
	my %ta = map { $_ => 1 } @fl_sizes; 
	( keys %ta ) == @fl_sizes or &stopErr("[Err] there are repeated flankSize.\n"); 
}
for (my $i=0; $i<@fl_sizes; $i++) {
	if (defined $opts{'flankName'}[$i] and $opts{'flankName'}[$i] ne '') {
		$fl_names[$i] = $opts{'flankName'}[$i]; 
	} else {
		$fl_names[$i] = $fl_sizes[$i]; 
	}
}
{
	my %ta = map { $_ => 1 } @fl_names; 
	(keys %ta) == @fl_names or &stopErr("[Err] there are repeated flankName.\n"); 
}
my $pl_dg = $opts{'pl_dg'} // '/home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl'; 
my $call_dg = ( -e $pl_dg ) ? "perl $pl_dg" : "$pl_dg" ; 

&tsmsg("[Rec] All begin.\n"); 
# Make non-overlap gene model files. 
my @out_files; 
my $nonOvl_file = 'non_ovl.combined'; 
for (my $i=0; $i<@bgGffs; $i++) {
	my $bg_gff = $bgGffs[$i]; 
	my $out_file = "non_ovl.$i"; 
	&exeCmd_1cmd("$call_dg -inGff $bg_gff -compare2gffC $opts{'seedGff'} -rmOvlap -rmOvlapLen 1 -rmOvlapType 'CDS,match_part' -rmOvlapStrand Single -out $out_file"); 
	push(@out_files, $out_file); 
}# End for 
&exeCmd_1cmd("cat $opts{'seedGff'} @out_files > $nonOvl_file"); 

# The previous non-overlap gene model file is $nonOvl_file. 
for (my $i=0; $i<@fl_sizes; $i++) {
	my $fl_len = $fl_sizes[$i]; 
	my $fl_tag = $fl_names[$i]; 
	my $island_file = "island.$i"; 
	my $good_file = "good.$i.$fl_tag"; 
	&exeCmd_1cmd("$call_dg -inGff $nonOvl_file -islandGene $fl_len -out $island_file -islandFeatType 'CDS,match_part' -islandStrand 'Both'"); 
	&exeCmd_1cmd("$call_dg -inGff $opts{seedGff} -compare2gffC $island_file -out $good_file -sameIntron -sameSingleExon"); 
	&exeCmd_1cmd("$call_dg -inGff $good_file -listTopID -out ${good_file}.ID"); 
}
&tsmsg("[Rec] All done.\n"); 

