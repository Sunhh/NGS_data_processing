#!/usr/bin/perl
# [3/17/2022] Need to add the cumulative coverage length of the query gene.
# [3/20/2022] Direction matters, so I should run bedtools twice with different directions.
# [3/22/2022] The output gff3 of liftoff needs to be filtered because there are still gene models with gene coverage lower than required.
#             Strand-specific overlapping has been fixed.
# [3/23/2022] I decide not to filter coverage in this program, but I'll output the coverage bases in the output intersection file.
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;

!@ARGV and die "perl $0 from_XX to_XX\n";

my $fromXX = shift;
my $toXX = shift;

my %info;
$info{'CL'} = ['ASM1002_97103'];
$info{'CM'} = ['ASM1004_USVL531'];
$info{'CA'} = ['ASM1005_PI482246'];
$info{'CC'} = ['ASM1003_PI537277'];
$info{'WCG'} = ['WCGv2'];
$info{'CORD'} = ['cordophanusV1.5'];

my %dta;
$dta{'para_minCov'}   = 0.5;
$dta{'para_minIdent'} = 0.75;
$dta{'para_liftoff'}  = "-p 20 -s $dta{'para_minIdent'} -a $dta{'para_minCov'} -copies";
$dta{'pl_cntR2Q'}   = 'perl /data/Sunhh/try_liftoff/tools/cnt_R2Q_liftoff_info.pl';
# $dta{'pl_cnvt_gff4igv'}  = 'perl /data/Sunhh/try_liftoff/tools//fit_gff_4igv.pl';
$dta{'exe_bedtools'}     = 'bedtools';
$dta{'exe_liftoff'}      = 'liftoff';
$dta{'tmp_dir'}  = "tmp_files/mapCDS.$info{$fromXX}[0].to.$info{$toXX}[0]/";
-e "tmp_files" or &runCmd("mkdir -p tmp_files/");
$dta{'chroms'}   = "input/map_chr_list/map.$info{$fromXX}[0].vs.$info{$toXX}[0]";
$dta{'unplaced'} = "input/map_chr_list/unplaced_seqID.$info{$fromXX}[0]";
$dta{'to_gff3'} = "input/$info{$toXX}[0].trim2CDS.gff3";
$dta{'from_gff3'} = "input/$info{$fromXX}[0].trim2CDS.gff3";
$dta{'from_fas'}  = "common_data/$info{$fromXX}[0].ref.fa";
$dta{'to_fas'}  = "common_data/$info{$toXX}[0].ref.fa";
$dta{'out_gff3'}     = "output/mapCDS.$info{$fromXX}[0].to.$info{$toXX}[0].gff3";
$dta{'out_tbl'}      = "output/mapCDS.$info{$fromXX}[0].to.$info{$toXX}[0].tbl";
$dta{'out_igvGff3'}  = "forIGV/mapCDS.$info{$fromXX}[0].to.$info{$toXX}[0].forIGV.gff3";
$dta{'out_unmap'}    = "output/mapCDS.$info{$fromXX}[0].to.$info{$toXX}[0].unmapped";
$dta{'toRM'} = [];

&runCmd("rm -rf $dta{'tmp_dir'}");
&runCmd("$dta{'exe_liftoff'} $dta{'para_liftoff'} -dir $dta{'tmp_dir'} -chroms $dta{'chroms'} -unplaced $dta{'unplaced'} -g $dta{'from_gff3'} -o $dta{'out_gff3'} -u $dta{'out_unmap'} $dta{'to_fas'} $dta{'from_fas'}");
&runCmd("$dta{'pl_cntR2Q'}   $dta{'from_gff3'}  $dta{'to_gff3'}  $dta{'out_gff3'} > $dta{'out_tbl'}");
push(@{$dta{'toRM'}}, $dta{'tmp_dir'});

for (@{$dta{'toRM'}}) {
  &fileSunhh::_rmtree($_);
}


