#!/usr/bin/perl
# 2019-03-08 Don't translate geno_tbl to struct. 
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"indv_txt:s", # col_[-3] is individual ID; 
	"snp_number:i", # 1000 
	"minR:i", # minimum K : 2
	"maxR:i", # maximum K : 7
	"geno_tbl:s", # No header line. 
	"isStructTbl!", 
	"OutDir:s", # pwd
	"BinDir:s", 
	"maxK:i", # 20 
	"minK:i", 
	"randSeed!", 
); 
my $curr_dir = `pwd`; chomp($curr_dir); 
$opts{'minR'} //= 1; 
$opts{'maxR'} //= $opts{'minR'} //= 20; 
$opts{'maxK'} //= 20; 
$opts{'minK'} //= 2; 
$opts{'snp_number'} //= 1000; 
$opts{'minR'} <= $opts{'maxR'} or die "minR > maxR\n"; 
$opts{'OutDir'} //= `pwd`; 
chomp($opts{'OutDir'}); 
$opts{'BinDir'} //= '/home/Sunhh/tools/github/NGS_data_processing/structure'; 
# die  "perl $0  <add_ref_list> <individual> <snp_number> <run>\n" unless ($#ARGV==2);
my $seed_tag = ''; 
$opts{'randSeed'} and $seed_tag = '-randSeed'; 

my $individual = $opts{'indv_txt'}; 
my $snp_number = $opts{'snp_number'}; 
my $UsePosition = $opts{'geno_tbl'}; 
my $UsePosition1 = $UsePosition; 
$UsePosition1 =~ s!^\S+/!!; 

my $OutDir = $opts{'OutDir'};

my $BinDir = $opts{'BinDir'}; 

sub usage {
	print STDERR <<HH; 

perl $0 -indv_txt indvID.shrtC-3 -snp_number 1000 -geno_tbl all.snp 

-minR / -maxR 
-minK / -maxK
-randSeed

-snp_number      [$opts{'snp_number'}] Number of SNP to used in each dataset. -1 means all of SNP will be used. 


HH
	exit(1); 
}

$opts{'help'} and &usage(); 
for (qw/indv_txt geno_tbl/) {
	defined $opts{$_} or do { print STDERR "Need -$_\n"; &usage(); }; 
}



for (my $i=$opts{'minR'};$i<=$opts{'maxR'};$i++){
	my  $out_Dir="$OutDir/structure_$i";
	mkdir($out_Dir) unless (-d $out_Dir); 
	chdir($out_Dir); 
	if ($snp_number > 0) {
		&exeCmd_1cmd("perl $BinDir/rand_small_position.pl $OutDir/$UsePosition  $out_Dir/$UsePosition1  $snp_number");
	} else {
		&exeCmd_1cmd("cp -p $OutDir/$UsePosition $out_Dir/$UsePosition1"); 
	}
	# my $location_num = `wc -l $out_Dir/$UsePosition1  | cut  -f 1 -d " " | tr -d '\n' `;
	my $c1 = "head -1 $out_Dir/$UsePosition1 | sed -e 's!\\s\\+!\\n!g' | tail -n +2 | wc -l | cut -f 1 -d \" \" | tr -d '\\n'"; 
	my $location_num = `$c1`; 
	my $sample_num   = `wc -l $curr_dir/$individual | cut  -f 1 -d " " | tr -d '\n'`; 
	

	if ($opts{'isStructTbl'}) {
		&exeCmd_1cmd("cp -p $out_Dir/$UsePosition1 $out_Dir/Input.file"); 
	} else {
		&exeCmd_1cmd("perl $BinDir/get_structure_input.pl $out_Dir/$UsePosition1 $curr_dir/$individual $out_Dir/Input.file");
	}
	&exeCmd_1cmd("perl $BinDir/new_mainparams.pl $seed_tag -minK $opts{'minK'} -K $opts{'maxK'} -output $out_Dir -input $BinDir  -loca $location_num  -sample $sample_num");

	# open O,'>', "cmd_list_runStruct" or die; 
	for (my $j=$opts{'minK'}; $j<=$opts{'maxK'}; $j++) {
		my $jj = sprintf("%02d", $j); 
		my $dd = "$out_Dir/structure_K${jj}"; 
		print STDOUT "cd $dd; $BinDir/structure 1>stdout 2>stderr ; cd -; \n"; 
	}
	# close O; 
	# `sh   $BinDir/Start.sh     $out_Dir  $i   `;
	chdir($OutDir); 
}
