#!/usr/bin/perl-w
use strict;
#explanation:this program is edited to 
#edit by HeWeiMing;   Wed Dec  8 18:03:40 HKT 2010
#Version 1.0    hewm@genomics.org.cn 

die  "Version 1.0\t2010-12-8;\nUsage: $0 <result.list><individual.txt><Out>\n" unless (@ARGV ==3);
#############Befor  Start  , open the files ####################
open (Res,"$ARGV[0]") || die "input file can't open $!";
open (IB,"<$ARGV[1]") || die "input file can't open $!";
open (OA,">$ARGV[2]") || die "output file can't open $!" ;
my %OUT ;
################ Do what you want to do #######################
while($_=<Res>)
{
    chomp ;
    open (IA,"<$_") or die "$!" ;
$/="Inferred clusters";
#Label (%Miss) :  Inferred clusters
<IA> ;
$/="\n";
$_=<IA>;
	while($_=<IA>) 
	{ 
        # 9     T_09    (3)   :  0.057 0.943
        #10     T_10    (0)   :  0.049 0.951
        #11     T_11    (0)   :  0.022 0.978
        chomp ; 
        last if  ($_ eq "" ) ;
		my @inf=split(/\:/ ,$_ );
        $inf[0]=~s/^\s+//; 
        $inf[-1]=~s/^\s+//; 
        my @aa=split /\s+/,$inf[0];
        my @bb=split /\s+/,$inf[-1];
        my $K=$#bb+1;
        my $Population=join("\t",@bb);
        $OUT{$aa[1]}=$OUT{$aa[1]}."\tK=$K\t$aa[1]\t$Population";
	}
close IA;
}
close Res ;
while(<IB>)
    {
        chomp  ;
        my @inf=split ;
        next if  (!exists $OUT{$inf[-3]} );
        print OA  "$inf[0]",$OUT{$inf[-3]},"\n";
    }
close OA ;
close IB ;
######################swimming in the sky and flying in the sea ###########################
