#!/usr/bin/perl-w
use strict;
#explanation:this program is edited to 
#edit by HeWeiMing;   Mon Sep 13 11:30:24 CST 2010
#Version 1.0    hewm@genomics.org.cn 

die  "Version 1.0\t2010-09-13;\nUsage: $0 <InPut_1><out><rand.cout>\n" unless (@ARGV >=2);

#############Befor  Start  , open the files ####################

open IA,"$ARGV[0]"  || die "input file can't open $!" ;

open OA,">$ARGV[1]" || die "output file can't open $!" ;

################ Do what you want to do #######################
my $cout=1; 
my @arry=();

while(<IA>) 
	{ 
     $arry[$cout]=$_ ;
     $cout++;    
#		chomp ; 
#		my @inf=split ;
	}
close IA;
    
    $ARGV[2]||=500000;
    if ($cout<$ARGV[2] )
    {
        print "bad"; exit(1);
    }
    my %hash=();
    for (my $ii=1; $ii<=$ARGV[2] ;  )
    {
        my $k=int(rand($cout))+1;
        if (!exists $hash{$k})
        {
            print OA  $arry[$k];
             $hash{$k}=1 ;
             $ii++;
        }
    }
close OA ;
    `sort  -s -k 1,1 -k 2n -T./  -o $ARGV[1]  $ARGV[1] ` ;
#    system ("mv $ARGV[1]  $ARGV[0] ") ;

#   sort -s -k 1,1 -k 2,2 -T./

######################swimming in the sky and flying in the sea ###########################
