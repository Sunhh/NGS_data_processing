#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"snp_tbl:s", # Acc131_mask.snp 
	"opref:s",   # 
	"tax_list:s",# taxa list 
	"replaceID:s",# We should change IDs in fasta 
	"wantMega!", # 
	"dirID:s",   # '01' 
	"setID:s",   # '01'
	"windN:i",   # 1
	"windL:i",   # 10000
	"pl_randSlct:s", # 'perl /home/Sunhh/tools/github/NGS_data_processing/reseq_tools/rand_site_wiWind.pl'
	"pl_dealTbl:s", 
	"pl_tbl2fas:s",  # 'perl /home/Sunhh/tools/github/NGS_data_processing/reseq_tools/cnvt_tools/tbl2fas.pl'
	"pl_dealFas:s",  # 'deal_fasta.pl'
	"pl_fas2meg:s",  # 'perl /home/Sunhh/tools/github/NGS_data_processing/reseq_tools/cnvt_tools/fas2meg.pl' 
); 

$opts{'windN'} //= 1; 
$opts{'windL'} //= 10000; 
$opts{'dirID'} //= '01'; 
$opts{'setID'} //= '01'; 
$opts{'pl_randSlct'} //= 'perl /home/Sunhh/tools/github/NGS_data_processing/reseq_tools/rand_site_wiWind.pl'; 
$opts{'pl_dealTbl'}  //= 'deal_table.pl'; 
$opts{'pl_tbl2fas'}  //= 'perl /home/Sunhh/tools/github/NGS_data_processing/reseq_tools/cnvt_tools/tbl2fas.pl'; 
$opts{'pl_dealFas'}  //= 'deal_fasta.pl'; 
$opts{'pl_fas2meg'}  //= 'perl /home/Sunhh/tools/github/NGS_data_processing/reseq_tools/cnvt_tools/fas2meg.pl'; 

my $help_txt = <<HH; 

perl $0 -snp_tbl Acc131_mask.snp   -opref Acc131_mask   -tax_list GrpList/grp52_list   -dirID '01' -setID '01' 

-help

-tax_list     [filename] First column is the header of snp table. 


-windN        [$opts{'windN'}]
-windL        [$opts{'windL'}]

-pl_randSlct  [$opts{'pl_randSlct'}]
-pl_dealTbl   [$opts{'pl_dealTbl'}]
-pl_tbl2fas   [$opts{'pl_tbl2fas'}]
-pl_dealFas   [$opts{'pl_dealFas'}]
-pl_fas2meg   [$opts{'pl_fas2meg'}]

HH

defined $opts{'snp_tbl'} or &LogInforSunhh::usage($help_txt); 
defined $opts{'opref'} or &LogInforSunhh::usage($help_txt); 
defined $opts{'tax_list'} or &LogInforSunhh::usage($help_txt); 

# # For set05 - 06_rand_05_30in10k/
# perl rand_site_wiWind.pl -snp_tbl Acc131_mask.snp -wind_len 10000 -wind_num 30 > 06_rand_05_30in10k/Acc131_mask_set05.snp
# deal_table.pl 06_rand_05_30in10k/Acc131_mask_set05.snp -colByTbl GrpList/grp52_list -colByTbl_also 0,1 > 06_rand_05_30in10k/Acc131_mask_set05.grp52.snp
# perl /home/Sunhh/tools/github/NGS_data_processing/reseq_tools/cnvt_tools/tbl2fas.pl -startColN 2   06_rand_05_30in10k/Acc131_mask_set05.grp52.snp | deal_fasta.pl -replaceID -replaceIDlist GrpList/taxaList_name -replaceIDcol 0,5 > 06_rand_05_30in10k/Acc131_mask_set05.grp52.snp.fa

my $setID = $opts{'setID'}; 
my $dirID = $opts{'dirID'}; 
my $wNum  = $opts{'windN'}; 
my $wLen  = $opts{'windL'}; 
my $opref = $opts{'opref'}; 

my $oDir = "${dirID}_rand_${setID}_${wNum}in" . int($wLen/1000) . "k"; 

&tsmsg("[Rec] Begin [$0]\n"); 
-d $oDir or mkdir($oDir); 
&exeCmd_1cmd( "$opts{'pl_randSlct'} -snp_tbl $opts{'snp_tbl'} -wind_len $wLen -wind_num $wNum > $oDir/${opref}_set${setID}.snp" ) and &stopErr("[Err]\n"); 
&exeCmd_1cmd( "$opts{'pl_dealTbl'} -colByTbl $opts{'tax_list'} -colByTbl_also 0,1 $oDir/${opref}_set${setID}.snp > $oDir/${opref}_set${setID}.use.snp" ) and &stopErr("[Err]\n"); 
&exeCmd_1cmd( "$opts{'pl_tbl2fas'} -startColN 2   $oDir/${opref}_set${setID}.use.snp > $oDir/${opref}_set${setID}.use.snp.ori.fa -showTime 100000" ) and &stopErr("[Err]\n"); 
if ( defined $opts{'replaceID'} ) {
	&exeCmd_1cmd( "$opts{'pl_dealFas'} -replaceID -replaceIDlist $opts{'tax_list'} -replaceIDcol $opts{'replaceID'} $oDir/${opref}_set${setID}.use.snp.ori.fa > $oDir/${opref}_set${setID}.use.snp.fa" ) and &stopErr("[Err]\n"); 
} else {
	&fileSunhh::_move( "$oDir/${opref}_set${setID}.use.snp.ori.fa", "$oDir/${opref}_set${setID}.use.snp.fa" ); 
}

if ($opts{'wantMega'}) {
	&exeCmd_1cmd("$opts{'pl_fas2meg'} $oDir/${opref}_set${setID}.use.snp.fa > $oDir/${opref}_set${setID}.use.snp.meg"); 
}

&tsmsg("[Rec] All done [$0]\n"); 

