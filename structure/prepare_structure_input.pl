#/usr/bin/perl
# 2019-03-08 Transform SNP tables into STRUCTURE input file.
use strict;
use warnings;
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"task:s",  
	  # cnvt_SNP 
	  # cnvt_vcfTab 
	  # get_indvList 
	  # to be added : shorten_indvList / get_indvList
	  # get_param
	"inFmt:s", # SNP / vcfTab
	"infile_SNP:s", # 
	"infile_indv:s", 
	"outfile_struct:s", 
	"outfile_sampleN:s", 
	"tmpDir:s", 
	# For shorten_indvList
	"shrt_len:i", # 9
	"shrt_col:i", # 0
	# For get_param_1
	### Required: 
	##### $opts{'outfile_struct'} or $opts{'infile_SNP'}
	##### $opts{'outfile_mainparam'}
	##### $opts{'outfile_extraparam'}
	"outfile_mainparam:s", 
	"outfile_extraparam:s", 
	"basic_mainparam:s", # /home/Sunhh/tools/github/NGS_data_processing/structure/mainparams
	"basic_extraparam:s", # /home/Sunhh/tools/github/NGS_data_processing/structure/extraparams
	"param_K:i",        # KKKKKKK; MAXPOPS 4
	"param_burn:i",     # AAAAAAA; BURNIN  40000
	"param_numrep:i",   # BBBBBBB; NUMREPS 15000

	"param_infile:s",   # CCCCCCC; INFILE  $opts{'outfile_struct'}
	"param_outfile:s",  # DDDDDDD; OUTFILE $opts{'outfile_struct'}.result

	"param_sampleN:i",  # EEEEEEE; NUMINDS From $opts{'outfile_struct'}
	"param_lociN:i",    # FFFFFFF; NUMLOCI From $opts{'outfile_struct'}

	"param_seedN:i",    # RRRRRRR; SEED    -1 means randomly select it. 
	"param_exe:s",      # /home/Sunhh/tools/github/NGS_data_processing/structure/structure
	# For get_param_n
	"param_minK:i",     # 1
	"param_maxK:i",     # 20
	"param_outDir:s",   # required. 
	"param_cmdLis:s",   # 
); 
my %globV; 
&set_glob(); 
$opts{'task'} eq 'CNVT_SNP' and &run_cnvt_SNP(); 
$opts{'task'} eq 'CNVT_VCFTAB' and &run_cnvt_SNP(); 
$opts{'task'} eq 'GET_INDVLIST' and &run_get_indvList(); 
$opts{'task'} eq 'GET_PARAM_1' and &run_get_param_1(); 
$opts{'task'} eq 'GET_PARAM_N' and &run_get_param_n(); 
sub run_get_param_n {
	defined $opts{'param_outDir'} or &LogInforSunhh::usage($globV{'htxt'}); 
	$opts{'param_outDir'} = &fileSunhh::_abs_path($opts{'param_outDir'}); 
	-e $opts{'param_outDir'} or mkdir($opts{'param_outDir'}); 
	-e $opts{'param_outDir'} or &stopErr("[Err] Failed to get dir [$opts{'param_outDir'}]\n"); 
	$opts{'param_cmdLis'} //= "$opts{'param_outDir'}/cmd_list"; 
	&fileSunhh::write2file($opts{'param_cmdLis'}, '', '>'); 
	$opts{'param_minK'} //= 1; 
	$opts{'param_maxK'} //= 20; 
	for (my $k=$opts{'param_minK'}; $k<=$opts{'param_maxK'}; $k++) {
		my $stru_dir = sprintf("%s/structure_K%02d", $opts{'param_outDir'}, $k); 
		$stru_dir = &fileSunhh::_abs_path($stru_dir); 
		-e $stru_dir or mkdir($stru_dir); 
		$opts{'param_K'} = $k; 
		$opts{'outfile_mainparam'}  = "$stru_dir/mainparams"; 
		$opts{'outfile_extraparam'} = "$stru_dir/extraparams"; 
		$opts{'param_outfile'}      = "$stru_dir/struct.result"; 
		my %ori; 
		for my $tk (qw/param_seedN/) {
			$ori{$tk} = $opts{$tk}; 
		}
		my $cmd = &run_get_param_1(); 
		for my $tk (qw/param_seedN/) {
			$opts{$tk} = $ori{$tk}; 
		}
		&fileSunhh::write2file($opts{'param_cmdLis'}, "cd $stru_dir/ ; $cmd 1>$stru_dir/stdout 2>$stru_dir/stderr\n", '>>'); 
	}

	return; 
}# run_get_param_n() 

defined $opts{'tmpDir'} and &fileSunhh::_rmtree($opts{'tmpDir'}); 
#### Subs 
sub run_get_param_1 {
	defined $opts{'basic_mainparam'} or $opts{'basic_mainparam'} = &_get_basic_mainparam(); 
	defined $opts{'basic_extraparam'} or $opts{'basic_extraparam'} = &_get_basic_extraparam(); 
	$opts{'param_exe'}    //= '/home/Sunhh/tools/github/NGS_data_processing/structure/structure'; 
	$opts{'param_K'}      //= 4; 
	$opts{'param_burn'}   //= 40000; 
	$opts{'param_numrep'} //= 15000; 
	$opts{'param_seedN'}  //= -1; 
	&_get_param_infile(); 
	&_get_param_outfile(); 
	&_get_param_sampleN(); 
	&_get_param_lociN(); 

	my $ifh_main = &openFH("$opts{'basic_mainparam'}", '<'); 
	my $ofh_main = &openFH("$opts{'outfile_mainparam'}", '>'); 
	while (<$ifh_main>) {
		s!KKKKKKK!$opts{'param_K'}!; 
		s!AAAAAAA!$opts{'param_burn'}!; 
		s!BBBBBBB!$opts{'param_numrep'}!; 
		s!CCCCCCC!$opts{'param_infile'}!; 
		s!DDDDDDD!$opts{'param_outfile'}!; 
		s!EEEEEEE!$opts{'param_sampleN'}!; 
		s!FFFFFFF!$opts{'param_lociN'}!; 
		print {$ofh_main} $_; 
	}
	close($ofh_main); 
	close($ifh_main); 

	my $ifh_extr = &openFH("$opts{'basic_extraparam'}", '<'); 
	my $ofh_extr = &openFH("$opts{'outfile_extraparam'}", '>'); 
	while (<$ifh_extr>) {
		if (m!RRRRRRR!) {
			my $seedN = $opts{'param_seedN'}; 
			$seedN > 0 or $seedN = &rand_9(); 
			s!RRRRRRR!$seedN!; 
		}
		print {$ofh_extr} $_; 
	}
	close($ofh_extr); 
	close($ifh_extr); 
	my $cmd = "$opts{'param_exe'} -m $opts{'outfile_mainparam'} -e $opts{'outfile_extraparam'}"; 
	&tsmsg("[Msg] to run cmd: $cmd 1>stdout 2>stderr\n"); 
	return($cmd); 
}# run_get_param_1() 

sub rand_9 {
	my $nn = ""; 
	for (1 .. 9) {
		my $r = int(rand(10)); 
		$r == 10 and $r = 9; 
		$nn .= $r; 
	}
	$nn =~ s!^0!1!; 
	return($nn); 
}# rand_9() 

sub _get_param_lociN {
	if (defined $opts{'param_lociN'} and $opts{'param_lociN'} > 0) {
		return; 
	}
	&_get_param_infile(); 
	# This is only good for two-row format. 
	my $ifh = &openFH($opts{'param_infile'}, '<'); 
	while (<$ifh>) {
		chomp; 
		s!^\s+!!; 
		my @ta = split(/\s+/, $_); 
		$opts{'param_lociN'} = scalar(@ta); 
		last; 
		# m!^(\S+)! and last; 
	}
	close($ifh); 
	return($opts{'param_lociN'}); 
}# _get_param_lociN() 
sub _get_param_sampleN {
	if (defined $opts{'param_sampleN'} and $opts{'param_sampleN'} > 0) {
		return; 
	}
	&_get_param_infile(); 
	my $ifh = &openFH($opts{'param_infile'}, '<'); 
	my $cnt = 0; 
	my %h; 
	while (<$ifh>) {
		m!^(\S+)! or next; 
		defined $h{$1} and next; 
		$h{$1} ++; 
		$cnt ++; 
	}
	close($ifh); 
	$opts{'param_sampleN'} = $cnt; 
	return($cnt); 
}# _get_param_sampleN() 
sub _get_param_outfile {
	if (defined $opts{'param_outfile'}) {
		$opts{'param_outfile'} = &fileSunhh::_abs_path($opts{'param_outfile'}); 
	} else {
		&_get_param_infile(); 
		$opts{'param_outfile'} = "$opts{'param_infile'}.result"; 
	}
	return; 
}# _get_param_outfile() 
sub _get_param_infile {
	if (defined $opts{'param_infile'}) {
		$opts{'param_infile'} = &fileSunhh::_abs_path($opts{'param_infile'}); 
	} elsif (defined $opts{'outfile_struct'}) {
		$opts{'param_infile'} //= &fileSunhh::_abs_path($opts{'outfile_struct'}); 
	} else {
		$opts{'param_infile'} //= 'input.struct'; 
	}
	return; 
}# _get_param_infile() 

sub run_get_indvList {
	$opts{'infile_indv'} //= "$opts{'infile_SNP'}.struct.indv"; 
	my $ifh_1 = &openFH("$opts{'infile_SNP'}",'<'); 
	my $h1 = <$ifh_1>; 
	chomp($h1); 
	my @h2 = &splitL("\t", $h1); 
	$opts{'inFmt'} eq 'VCFTAB' and splice(@h2, 2, 1); 
	splice(@h2, 0, 2); 
	close($ifh_1); 
	my $ofh_1 = &openFH("$opts{'infile_indv'}", '>'); 
	my $fakeV = 0; 
	my %h; 
	for my $a1 (@h2) {
		$fakeV ++; 
		my $fakeID = sprintf("I%04d", $fakeV); 
		my $rawID = $a1; 
		defined $h{$fakeID} and die "Failed to find fakeID for [$a1]\n"; 
		$h{$fakeID} = 1; 
		print {$ofh_1} join("\t", $fakeID, $rawID, $rawID)."\n"; 
	}
	close($ofh_1); 
	return; 
}# run_get_indvList() 

sub run_cnvt_SNP {
	defined $opts{'infile_SNP'} or defined $opts{'infile_indv'} or &LogInforSunhh::usage($globV{'htxt'}); 
	$opts{'outfile_struct'} //= "$opts{'infile_SNP'}.struct"; 
	$opts{'outfile_sampleN'} //= "$opts{'outfile_struct'}.sampleN"; 
	$opts{'infile_indv'}    //= "$opts{'outfile_struct'}.indv"; 
	-e $opts{'infile_indv'} or &run_get_indvList(); 
	
	my @name_list = map { $_->[0] } &fileSunhh::load_tabFile($opts{'infile_indv'}); 
	my @out_rows; 
	push(@out_rows, [' ']); 
	push(@out_rows, map { [$_, $_] } @name_list); 
	$globV{'sample_num'} = scalar(@name_list); 
	
	my $ifh_snp = &openFH($opts{'infile_SNP'}, '<'); 
	my ($prevChr, $prevPos) = ('', ''); 
	while (<$ifh_snp>) {
		chomp; 
		my @ta = &splitL("\t", $_); 
		$opts{'inFmt'} eq 'VCFTAB' and splice(@ta, 2, 1); 
		if ($. == 1 and $ta[0] =~ m!^\#?(chr|chrom|chromosome)$!i) {
			scalar(@ta)-2 == $globV{'sample_num'} or &stopErr("[Err] Bad sample number [$globV{'sample_num'}]: @ta\n"); 
			next; 
		}
		if ($prevChr ne $ta[0]) {
			$out_rows[0][0] .= " -1"; 
			$prevChr = $ta[0]; 
			$prevPos = $ta[1]; 
		} else {
			$out_rows[0][0] .= (" " . ($ta[1]-$prevPos)); 
		}

		my %uniqB2V; $uniqB2V{'maxV'} = 100; 
		if ($opts{'inFmt'} eq 'VCFTAB') {
			for (my $i=2; $i<@ta; $i++) {
				$ta[$i] =~ m!^(\S+)/(\S+)$! or &stopErr("[Err] Bad genotype [$ta[$i]]\n"); 
				my @aa = ($1, $2); 
				for (my $j=0; $j<@aa; $j++) {
					if (defined $globV{'iupac'}{$aa[$j]}) {
						$out_rows[$i-1][$j] .= " $globV{'iupac'}{$aa[$j]}[0]"; 
					} elsif (defined $uniqB2V{'allele'}{$aa[$j]}) {
						$out_rows[$i-1][$j] .= " $uniqB2V{'allele'}{$aa[$j]}[0]"; 
					} else {
						$uniqB2V{'maxV'} ++; 
						$uniqB2V{'allele'}{$aa[$j]} = [$uniqB2V{'maxV'}, $uniqB2V{'maxV'}]; 
						&tsmsg("[Wrn] Assigned a new number to allele [$aa[$j]]\n"); 
						$out_rows[$i-1][$j] .= " $uniqB2V{'allele'}{$aa[$j]}[0]"; 
					}
				}
			}
		} elsif ($opts{'inFmt'} eq 'SNP') {
			for (my $i=2; $i<@ta; $i++) {
				if (defined $globV{'iupac'}{$ta[$i]}) {
					$out_rows[$i-1][0] .= " $globV{'iupac'}{$ta[$i]}[0]"; 
					$out_rows[$i-1][1] .= " $globV{'iupac'}{$ta[$i]}[1]"; 
				} elsif (defined $uniqB2V{'allele'}{$ta[$i]}) {
					$out_rows[$i-1][0] .= " $uniqB2V{'allele'}{$ta[$i]}[0]"; 
					$out_rows[$i-1][1] .= " $uniqB2V{'allele'}{$ta[$i]}[1]"; 
				} else {
					$uniqB2V{'maxV'} ++; 
					$uniqB2V{'allele'}{$ta[$i]} = [$uniqB2V{'maxV'}, $uniqB2V{'maxV'}]; 
					&tsmsg("[Wrn] Assigned a new number to allele [$ta[$i]]\n"); 
					$out_rows[$i-1][0] .= " $uniqB2V{'allele'}{$ta[$i]}[0]"; 
					$out_rows[$i-1][1] .= " $uniqB2V{'allele'}{$ta[$i]}[1]"; 
				}
			}
		}
	}
	close($ifh_snp); 
	&fileSunhh::write2file( $opts{'outfile_struct'}, join("\n", map { @$_ } @out_rows)."\n", '>' ); 
	&fileSunhh::write2file( $opts{'outfile_sampleN'}, "$globV{'sample_num'}\n", '>'); 
}# run_cnvt_SNP() 


sub set_glob {
	$globV{'htxt'} = <<HH; 
################################################################################
# perl $0     -task [cnvt_SNP/cnvt_vcfTab]   -inFmt [SNP/vcfTab]   -infile_SNP in.snp   -outfile_struct in.snp.struct
################################################################################
HH
	$opts{'help'} and &LogInforSunhh::usage($globV{'htxt'}); 
	$opts{'task'}   //= "cnvt_SNP"; 
	$opts{'inFmt'}  //= "SNP"; 
	$globV{'iupac'} = {
		'M' => [1,3], 
		'K' => [4,2], 
		'Y' => [3,2], 
		'R' => [1,4], 
		'W' => [1,2], 
		'S' => [3,4], 
		'A' => [1,1], 
		'T' => [2,2], 
		'C' => [3,3], 
		'G' => [4,4], 
		'-' => [-9,-9], 
		'.' => [-9,-9], 
		'N' => [-9,-9], 
	}; 
	$opts{'task'}  = uc($opts{'task'}); 
	$opts{'inFmt'} = uc($opts{'inFmt'}); 
	if ($opts{'task'} eq 'CNVT_VCFTAB') {
		$opts{'inFmt'} = 'VCFTAB'; 
	}
	# $opts{'tmpDir'} //= &fileSunhh::new_tmp_dir('create' => 1); 
}# 

sub _get_basic_extraparam {
	defined $opts{'tmpDir'} or $opts{'tmpDir'} = &fileSunhh::new_tmp_dir('create'=>1); 
	my $txt = <<'EPARAM'; 

EXTRA PARAMS FOR THE PROGRAM structure.  THESE PARAMETERS CONTROL HOW THE
PROGRAM RUNS.  ATTRIBUTES OF THE DATAFILE AS WELL AS K AND RUNLENGTH ARE 
SPECIFIED IN mainparams.

"(int)" means that this takes an integer value.
"(d)"   means that this is a double (ie, a Real number such as 3.14).
"(B)"   means that this variable is Boolean 
        (ie insert 1 for True, and 0 for False).

PROGRAM OPTIONS

#define NOADMIX     0 // (B) Use no admixture model (0=admixture model, 1=no-admix)
#define LINKAGE     0 // (B) Use the linkage model model 
#define USEPOPINFO  0 // (B) Use prior population information to pre-assign individuals
                             to clusters
#define LOCPRIOR    0 //(B)  Use location information to improve weak data

#define FREQSCORR   1 // (B) allele frequencies are correlated among pops
#define ONEFST      0 // (B) assume same value of Fst for all subpopulations.

#define INFERALPHA  1 // (B) Infer ALPHA (the admixture parameter)
#define POPALPHAS   0 // (B) Individual alpha for each population
#define ALPHA     1.0 // (d) Dirichlet parameter for degree of admixture 
                             (this is the initial value if INFERALPHA==1).

#define INFERLAMBDA 0 // (B) Infer LAMBDA (the allele frequencies parameter)
#define POPSPECIFICLAMBDA 0 //(B) infer a separate lambda for each pop 
					(only if INFERLAMBDA=1).
#define LAMBDA    1.0 // (d) Dirichlet parameter for allele frequencies 




PRIORS

#define FPRIORMEAN 0.01 // (d) Prior mean and SD of Fst for pops.
#define FPRIORSD   0.05  // (d) The prior is a Gamma distribution with these parameters

#define UNIFPRIORALPHA 1 // (B) use a uniform prior for alpha;
                                otherwise gamma prior
#define ALPHAMAX     10.0 // (d) max value of alpha if uniform prior
#define ALPHAPRIORA   1.0 // (only if UNIFPRIORALPHA==0): alpha has a gamma 
                            prior with mean A*B, and 
#define ALPHAPRIORB   2.0 // variance A*B^2.  


#define LOG10RMIN     -4.0   //(d) Log10 of minimum allowed value of r under linkage model
#define LOG10RMAX      1.0   //(d) Log10 of maximum allowed value of r
#define LOG10RPROPSD   0.1   //(d) standard deviation of log r in update
#define LOG10RSTART   -2.0   //(d) initial value of log10 r

                         
USING PRIOR POPULATION INFO (USEPOPINFO)

#define GENSBACK    2  //(int) For use when inferring whether an indiv-
                         idual is an immigrant, or has an immigrant an-
                         cestor in the past GENSBACK generations.  eg, if 
                         GENSBACK==2, it tests for immigrant ancestry 
                         back to grandparents. 
#define MIGRPRIOR 0.01 //(d) prior prob that an individual is a migrant 
                             (used only when USEPOPINFO==1).  This should 
                             be small, eg 0.01 or 0.1.
#define PFROMPOPFLAGONLY 0 // (B) only use individuals with POPFLAG=1 to update	P.
                                  This is to enable use of a reference set of 
                                  individuals for clustering additional "test" 
                                  individuals.

LOCPRIOR MODEL FOR USING LOCATION INFORMATION

#define LOCISPOP      0    //(B) use POPDATA for location information 
#define LOCPRIORINIT  1.0  //(d) initial value for r, the location prior
#define MAXLOCPRIOR  20.0  //(d) max allowed value for r




OUTPUT OPTIONS

#define PRINTNET     1 // (B) Print the "net nucleotide distance" to screen during the run
#define PRINTLAMBDA  1 // (B) Print current value(s) of lambda to screen
#define PRINTQSUM    1 // (B) Print summary of current population membership to screen

#define SITEBYSITE   0  // (B) whether or not to print site by site results. 
		     	       (Linkage model only) This is a large file!
#define PRINTQHAT    0  // (B) Q-hat printed to a separate file.  Turn this 
                           on before using STRAT.
#define UPDATEFREQ   100  // (int) frequency of printing update on the screen.
                                 Set automatically if this is 0.
#define PRINTLIKES   0  // (B) print current likelihood to screen every rep
#define INTERMEDSAVE 0  // (int) number of saves to file during run

#define ECHODATA     1  // (B) Print some of data file to screen to check
                              that the data entry is correct.
(NEXT 3 ARE FOR COLLECTING DISTRIBUTION OF Q:)
#define ANCESTDIST   0  // (B) collect data about the distribution of an-
                              cestry coefficients (Q) for each individual
#define NUMBOXES   1000 // (int) the distribution of Q values is stored as 
                              a histogram with this number of boxes. 
#define ANCESTPINT 0.90 // (d) the size of the displayed probability  
                              interval on Q (values between 0.0--1.0)



MISCELLANEOUS

#define COMPUTEPROB 1     // (B) Estimate the probability of the Data under 
                             the model.  This is used when choosing the 
                             best number of subpopulations.
#define ADMBURNIN  500    // (int) [only relevant for linkage model]: 
                             Initial period of burnin with admixture model (see Readme)
#define ALPHAPROPSD 0.025 // (d) SD of proposal for updating alpha
#define STARTATPOPINFO 0  // Use given populations as the initial condition 
                             for population origins.  (Need POPDATA==1).  It 
                             is assumed that the PopData in the input file 
                             are between 1 and k where k<=MAXPOPS.
#define RANDOMIZE      0  // (B) use new random seed for each run 
#define SEED     RRRRRRR  // (int) seed value for random number generator 
	                     (must set RANDOMIZE=0) 
#define METROFREQ    10   // (int) Frequency of using Metropolis step to update
                             Q under admixture model (ie use the metr. move every
                             i steps).  If this is set to 0, it is never used.
                             (Proposal for each q^(i) sampled from prior.  The 
                             goal is to improve mixing for small alpha.)
#define REPORTHITRATE 0 //   (B) report hit rate if using METROFREQ
EPARAM
	&fileSunhh::write2file("$opts{'tmpDir'}/extraparams", $txt, '>'); 
	return("$opts{'tmpDir'}/extraparams"); 
}# _get_basic_extraparam() 

sub _get_basic_mainparam {
	defined $opts{'tmpDir'} or $opts{'tmpDir'} = &fileSunhh::new_tmp_dir('create'=>1); 
my $txt = <<'MPARAM'; 

KEY PARAMETERS FOR THE PROGRAM structure.  YOU WILL NEED TO SET THESE
IN ORDER TO RUN THE PROGRAM.  VARIOUS OPTIONS CAN BE ADJUSTED IN THE
FILE extraparams.


"(int)" means that this takes an integer value.
"(B)"   means that this variable is Boolean 
        (ie insert 1 for True, and 0 for False)
"(str)" means that this is a string (but not enclosed in quotes!) 


Basic Program Parameters

#define MAXPOPS   KKKKKKK // (int) number of populations assumed
#define BURNIN    AAAAAAA // (int) length of burnin period
#define NUMREPS   BBBBBBB // (int) number of MCMC reps after burnin

Input/Output files

#define INFILE   CCCCCCC  // (str) name of input data file
#define OUTFILE  DDDDDDD  //(str) name of output data file

Data file format

#define NUMINDS  EEEEEEE  // (int) number of diploid individuals in data file
#define NUMLOCI  FFFFFFF  // (int) number of loci in data file
#define PLOIDY       2    // (int) ploidy of data
#define MISSING     -9    // (int) value given to missing genotype data
#define ONEROWPERIND 0    // (B) store data for individuals in a single line


#define LABEL     1     // (B) Input file contains individual labels
#define POPDATA   0     // (B) Input file contains a population identifier
#define POPFLAG   0     // (B) Input file contains a flag which says 
                              whether to use popinfo when USEPOPINFO==1
#define LOCDATA   0     // (B) Input file contains a location identifier

#define PHENOTYPE 0     // (B) Input file contains phenotype information
#define EXTRACOLS 0     // (int) Number of additional columns of data 
                             before the genotype data start.

#define MARKERNAMES      0  // (B) data file contains row of marker names
#define RECESSIVEALLELES 0  // (B) data file contains dominant markers (eg AFLPs)
                            // and a row to indicate which alleles are recessive
#define MAPDISTANCES     1  // (B) data file contains row of map distances 
                            // between loci


Advanced data file options

#define PHASED           0 // (B) Data are in correct phase (relevant for linkage model only)
#define PHASEINFO        0 // (B) the data for each individual contains a line
                                  indicating phase (linkage model)
#define MARKOVPHASE      0 // (B) the phase info follows a Markov model.
#define NOTAMBIGUOUS  -999 // (int) for use in some analyses of polyploid data



Command line options:

-m mainparams
-e extraparams
-s stratparams
-K MAXPOPS 
-L NUMLOCI
-N NUMINDS
-i input file
-o output file
-D SEED

MPARAM
	&fileSunhh::write2file("$opts{'tmpDir'}/mainparams", $txt, '>'); 
	return("$opts{'tmpDir'}/mainparams"); 
}# _get_basic_mainparam() 



