#/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

!@ARGV and die "perl $0 wind_list grpList_1 grpList_2 work_dir\n"; 

my $oDir = pop(@ARGV); 
-e $oDir and die "exist $oDir\n"; 
mkdir($oDir); 

open F,'<',"$ARGV[0]" or die; 
my %tmp_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>5e5 ); 
while (<F>) {
	&fileSunhh::log_section( $. , \%tmp_cnt ) and &tsmsg( "[Msg] $. line.\n" ); 
	chomp; 
	my @ta=split(/\t/, $_); 
	$ta[1] =~ m/^chr(\d+)$/i or next; 
	my $cn = $1; 
	my $oWind = $ta[0]; 
	$oWind =~ m!^.*/([^/]+)$! or die "bad $oWind\n"; 
	$oWind = "$oDir/$1"; 
	&fileSunhh::_copy( $ta[0], $oWind ); 
	print ("perl prepare_xpclr_input_wiGmP.pl $oWind.xpclr $oWind $ARGV[1] $ARGV[2] ; XPCLR -xpclr $oWind.xpclr.chr${cn}_g1.geno $oWind.xpclr.chr${cn}_g2.geno $oWind.xpclr.chr${cn}.snp $oWind.xpclr.chr${cn}.snp.out -w1 0.0005 100 100 $cn -p0 0.7 ;\n"); 
}
close F; 

