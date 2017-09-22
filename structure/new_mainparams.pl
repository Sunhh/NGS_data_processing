#!/usr/bin/perl -w
use strict;
#use Data::Dumper;  
#hewm@genomics.org.cn
use Getopt::Long;

####################
##USAGE
####################

sub usage{
	print STDERR <<USAGE;
Version:1.0
8/19/2009       hewm\@genomics.org.cn

		Usage: $0  
		where: new the mainparams of perl Structure
		      
	    Options
		   -input   <c> : InPutDir, Structure Temple[the temple mainparams] 
		   -output  <c> : OutPutDir ,[default ../ ]
		   -sample  <n> : the number of the sample (Need)
		   -loca    <n> : the number of the snp (need)
		   -K       <n> : number of populations assumed ,max K[default 7 ]
		   -burn    <n> : length of burnin period [default 30000]
		   -filein  <c> : Structure input files [default Input.file ]
		   -fileout <c> : Structure output file  [default result# ]
		   -numrep  <n> : number of MCMC reps after burning [10000]
		   -randSeed    : Generate random seed number in extraparams. 
		   -h           : show this help message
		   
USAGE
}

my ($input ,$sample,$loca,$K, $minK,$burn,$numrep,$help,$output, $filein,$fileout, $randSeed) ;

GetOptions(
	"input:s"=>\$input,
	"sample:s"=>\$sample,
	"loca:s"=>\$loca,
	"K:s"=>\$K,
	"minK:i"=>\$minK, 
	"filein:s"=>\$filein,
	"fileout:s"=>\$fileout,
	"burn:s"=>\$burn,
	"numrep:s"=>\$numrep,
	"output:s"=>\$output,
	"help"=>\$help,
	"randSeed!"=>\$randSeed, 
);

##check parameters and write into memory

if(defined($help)  || !defined($sample) || !defined($loca))
{
	usage;
	exit;
}

$K  //= 20 ;
$minK //= 2; $minK >= 2 or $minK = 2; 
$numrep //= 15000 ;
$burn  //=40000;
$fileout //="result" ;
$filein //="Input.file";
my $seedNum = 2245; 

 my  $TempDir=`pwd`;
 chomp  $TempDir ; 
 $input ||= $TempDir ;
 my @aa=split /\//,$TempDir ; 
 $TempDir=join("\/",@aa[0..$#aa-1]);
 $output ||=$TempDir ;
 my $CP=$input."\/mainparams" ;
 my $CE=$input."\/extraparams"; 
 
#############new the  script  ###########
for(my $ii=$minK; $ii<=$K; $ii++)
{
	my $aa=$ii ; 
	if($aa<10) {$aa="0".$aa ;}
	my $outDir=$output."/structure_K$aa" ;
	
	# system ("cp -r $input/structure_Temple  $outDir ") unless (  -e $outDir  ) ; 
	-e $outDir or mkdir($outDir); 
	system ("cp $filein $outDir/"); 
	my $outmainparams=$outDir."/mainparams";
	my $outextrparams=$outDir."/extraparams";
	my $result="$fileout".$aa ;
	
	open IN , "$CP"  || die "Can't open the Temple $!";
	open OUT , ">$outmainparams"  || die "Can't Write the Temple Here $!";
	
	while(<IN>)
	{
		s/KKKKKKK/$ii/ ;
		s/AAAAAAA/$burn/ ;
		s/BBBBBBB/$numrep/;
		s/CCCCCCC/$filein/;
		s/DDDDDDD/$result/;
		s/EEEEEEE/$sample/;
		s/FFFFFFF/$loca/;
		print OUT $_;    
	}
	close OUT;
	close IN;
	open CE,'<',"$CE" or die "Failed CE[$CE]\n"; 
	open OE,'>',"$outextrparams" or die "$!\n"; 
	while (<CE>) {
		if (m/RRRRRRR/) {
			$randSeed and $seedNum = int(rand(1000000));
			s/RRRRRRR/$seedNum/;
		}
		print OE $_; 
	}
	close OE; 
	close CE; 

	print "\t\tstructure_K$aa\thave\tdone\n";
}

######## swimming in the sky and flying in the sea#####
