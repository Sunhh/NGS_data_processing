#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"ref_goslim:s", # goslim_plant.obo.20181129
	"ref_bgobo:s",  # gene_ontology_edit.obo.2018-05-01 
	"ref_bggaf:s",  # ngsV1.GAF 
	"in_idlist:s",  # ID list; 

	"opref:s",      # out_prefix

	"pl_obo2tab:s", # /home/Sunhh/tools/github/NGS_data_processing/enrich/cnvt_GOobo_to_tab.pl
	"pl_map2slim:s", # map2slim
); 

defined $opts{'in_idlist'} or die "perl $0 -in_idlist in_geneID.list   -opref out_prefix   -ref_goslim goslim_plant.obo.20181129  -ref_bgobo gene_ontology_edit.obo.2018-05-01  -ref_bggaf ngsV1.GAFngsV1.GAF\n"; 

$opts{'pl_obo2tab'} //= 'perl /home/Sunhh/tools/github/NGS_data_processing/enrich/cnvt_GOobo_to_tab.pl'; 
$opts{'pl_map2slim'} //= 'map2slim'; 
$opts{'ref_goslim'}  //= '/home/Sunhh/tools/github/NGS_data_processing/enrich/goslim_plant.obo.20181129'; 
# $opts{'ref_bgobo'}   //= '/home/Sunhh/tools/github/NGS_data_processing/enrich/gene_ontology_edit.obo.2018-05-01.gz'; 
$opts{'ref_bgobo'}   //= '/home/Sunhh/tools/github/NGS_data_processing/enrich/gene_ontology_edit.obo.2018-05-01.gz'; 
$opts{'ref_bggaf'}   //= '/Data/Sunhh/watermelon/rnaseq/01_work/ref_WM97ngsV1Scf/for_fruitDev/ref/ngsV1.GAF'; 
$opts{'opref'}       //= 'opref'; 

my $wdir = &fileSunhh::new_tmp_dir( 'create'=>1 ); 

if ($opts{'ref_bgobo'} =~ m!^\.gz$!) {
	&runCmd("gzip -cd $opts{'ref_bgobo'} > $wdir/ref_bgobo"); 
	$opts{'ref_bgobo'} = "$wdir/ref_bgobo"; 
}
&runCmd("$opts{'pl_obo2tab'} $opts{'ref_goslim'} > $wdir/ref_goslim.tab"); 
my %cls_slim; 
{
	my @v1 = &fileSunhh::load_tabFile("$wdir/ref_goslim.tab"); 
	for (my $i=0; $i<@v1; $i++) {
		my @v2 = @{$v1[$i]}; 
		$cls_slim{$v2[0]} = [ $i, @v2[2,3,4] ]; 
	}
}

my %wantID = map { $_->[0] => 1} &fileSunhh::load_tabFile($opts{'in_idlist'}); 
my $fh_i1 = &openFH($opts{'ref_bggaf'}, '<'); 
my $fh_o1 = &openFH("$wdir/slct_gaf", '>'); 
while (<$fh_i1>) {
	chomp; 
	m/^\s*\!/ and do { print {$fh_o1} "$_\n"; next; }; 
	m!^\s*$! and next; 
	my @ta=&splitL("\t", $_); 
	defined $wantID{$ta[1]} or next; 
	print {$fh_o1} "$_\n"; 
}
close ($fh_o1); 
close ($fh_i1); 

# map2slim goslim_plant.obo.20181129 gene_ontology_edit.obo.2018-05-01 tt.gaf > tt.gaf.slim
&runCmd("$opts{'pl_map2slim'} $opts{'ref_goslim'} $opts{'ref_bgobo'} $wdir/slct_gaf > $wdir/slct_gaf.slim"); 
open F2,'<',"$wdir/slct_gaf.slim" or die; 
my $ofh_2 = &openFH("$opts{'opref'}.gene2class", '>'); 
my $ofh_3 = &openFH("$opts{'opref'}.cntClass", '>'); 
my %cnt_cls; 
print {$ofh_2} join("\t", qw/Gene_ID GO_ID Type GO_name/)."\n"; 
while (<F2>) {
	m/^\s*\!/ and next; 
	chomp; 
	my @ta = &splitL("\t", $_); 
	defined $cls_slim{$ta[4]} or die "failed to find GO_ID [$ta[4]]\n"; 
	defined $cnt_cls{$ta[4]}{$ta[1]} and next; 
	print {$ofh_2} join("\t", $ta[1], $ta[4], @{$cls_slim{$ta[4]}}[1,2])."\n"; 
	$cnt_cls{$ta[4]}{$ta[1]} ++; 
}
close F2; 
print {$ofh_3} join("\t", qw/Category GO_ID Gene_Num GO_name GO_def Gene_List/)."\n"; 
for my $k1 (sort { $cls_slim{$a}[1] cmp $cls_slim{$b}[1] || $cls_slim{$a}[0] <=> $cls_slim{$b}[0] } keys %cnt_cls) {
	my @allG = sort keys %{$cnt_cls{$k1}}; 
	my $n = scalar(@allG); 
	print {$ofh_3} join("\t", 
		$cls_slim{$k1}[1], 
		$k1, 
		$n, 
		$cls_slim{$k1}[2], 
		$cls_slim{$k1}[3], 
		join(",", @allG)
	)."\n"; 
}
close($ofh_2); 
close($ofh_3); 

&fileSunhh::_rmtree($wdir); 

sub runCmd {
	&exeCmd_1cmd($_[0]) and &stopErr("[Err] Failed at cmd: $_[0]\n"); 
}
