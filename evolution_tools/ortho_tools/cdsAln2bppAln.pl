#!/usr/bin/perl
# 9/1/2022: Generate sequence_alignment.txt and control.txt files for BP&P (bpp) software.
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;
use fastaSunhh;
my $fas_obj=fastaSunhh->new();
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "in_imap:s",    # Required.
  "in_cdsList:s", # Required.
  "out_prefix:s", # Default: bppIn.
  "help!",
);

my $htxt = <<HH;
####################################################################################################
# perl $0  -in_imap  bpp/OG.Imap.txt  -in_cdsList bpp/alnCds_list  -out_prefix bpp/bppOut.
#
# 
# -help     [Boolean]
#
####################################################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($htxt);
for my $k1 (qw/in_imap in_cdsList/) {
  defined $opts{$k1} or &LogInforSunhh::usage($htxt);
  $opts{$k1} = &fileSunhh::_abs_path_4link($opts{$k1});
}

$opts{'out_prefix'} //= 'bppOut.';

my $ofn1 = &fileSunhh::_abs_path_4link($opts{'out_prefix'} . "aln.txt");
my $ofn2 = &fileSunhh::_abs_path_4link($opts{'out_prefix'} . "bpp.ctl");
&fileSunhh::write2file($ofn1, '', '>');
&fileSunhh::write2file($ofn2, '', '>');
my $bppOutFn  = &fileSunhh::_abs_path_4link($opts{'out_prefix'} . "out.txt");
my $bppMcmcFn = &fileSunhh::_abs_path_4link($opts{'out_prefix'} . "mcmc.txt");

my (@seq_IDs, @spec_IDs);
for my $k1 (&fileSunhh::load_tabFile($opts{'in_imap'})) {
  push(@seq_IDs,  $k1->[0]);
  push(@spec_IDs, $k1->[1]);
}
my $seq_num = scalar(@seq_IDs);
my $spec_tree_num_line = join(" ", ('1') x $seq_num);
my $spec_tree_phase    = join(" ", ('0') x $seq_num);
my $spaceLen = 20;
for my $t1 (@seq_IDs) {
  $spaceLen < length($t1)+2 and $spaceLen = length($t1)+2;
}

my $aln_num = 0;
for my $l1 (&fileSunhh::load_tabFile($opts{'in_cdsList'})) {
  my ($fn1, $fnTag) = @$l1;
  $fnTag //= "";
  my %s = %{$fas_obj->save_seq_to_hash('faFile' => $fn1)};
  $aln_num ++;
  my $seq_len;
  for my $k1 (keys %s) {
    $s{$k1}{'seq'} =~ s!\s!!g;
    $seq_len //= length($s{$k1}{'seq'});
    $seq_len == length($s{$k1}{'seq'}) or &stopErr("[Err] sequence length varies in alignment [$fn1]\n");
    while ($s{$k1}{'seq'} =~ s!^(\?*)\-!$1?!) {;}
    while ($s{$k1}{'seq'} =~ s!\-(\?*)$!?$1!) {;}
    # $s{$k1}{'seq'} =~ s!\-!?!g;
  }
  &fileSunhh::write2file($ofn1, "$seq_num $seq_len\n\n", '>>');
  for my $taxID (@seq_IDs) {
    defined $s{$taxID} or &stopErr("[Err] Failed to find sequence for tax ID [$taxID]\n");
    &fileSunhh::write2file($ofn1, $fnTag.sprintf("^%-${spaceLen}s%s\n", $taxID, $s{$taxID}{'seq'}), '>>');
  }
  &fileSunhh::write2file($ofn1, "\n", '>>');
}

my $imap_txt= <<IMAP;
          seed =  -1

       seqfile = $ofn1
      Imapfile = $opts{'in_imap'}
       outfile = $bppOutFn
      mcmcfile = $bppMcmcFn

  speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2    * species delimitation rjMCMC algorithm0 and finetune(e)
* speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)

          speciestree = 0 * species tree NNI/SPR
*         speciestree = 1  0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio
*         speciestree = 1  * NNI over species/guide trees
*   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

         species&tree = $seq_num @spec_IDs
                           $spec_tree_num_line
                        ((K, C), (L, H));
                phase = $spec_tree_phase * 0: (default) indicating that sequences from that species are fully phased 
                                         *    haplotype sequences, so do not need to phase.
                                         * 1: unphased diploid sequences, so need to phase diploid unphased sequences.
   
              usedata = 1  * 0: no data (prior); 1:seq like
                nloci = $aln_num  * number of data sets in seqfile
              * nloci = 50  * Just a try run.

            cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?

           thetaprior = gamma 2 2000 # gamma(a, b) for theta (estimate theta)
         * thetaprior = 2 8.66  # gamma(a, b) for theta
                                #   Little is known for theta value.
                                #   Identity% between Walnut and others is about 76.91%;
                                #   So mean_theta=1-0.7691=0.2309;
                                #   Make a=2, then b=a/m=8.661758;
                                # [9/1/2022] This should the parameter set used in Maxima and Moschata genome paper.
             tauprior = gamma 2 1000 # gamma(a, b) for root tau & Dirichlet(a) for other tau's
           * tauprior = 282 2445  # gamma(a, b) for root tau & Dirichlet(a) for other tau's
                                  #   Identity% between Walnut and others is about 76.91%; Divergence time ~ 84 Mya;
                                  #   Then the mutation rate is about 0.2309/(84e6) /2 =~ 1.374405e-09 ;
                                  #   So the mean_tau(wn_other) = m = 84e6 * 1.374405e-09 = 0.11545 ;
                                  #   If make time +- 10Myr, then s=5e6 * 1.374405e-09 = 0.006872025; a=(m/s)**2=~282.2399, b=(m/s**2)=~2444.694;

      *      heredity = 0 4 4   # (0: No variation, 1: estimate, 2: from file) & a_gamma b_gamma (if 1)
      *     locusrate = 0 2.0   # (0: No variation, 1: estimate, 2: from file) & a_Dirichlet (if 1)
      * sequenceerror = 0 0 0 0 : 0.05 1   # sequencing errors: a_gamma, b_gamma. model of sequencing errors has changed, to be described later

      *      finetune = 0: 0.5 0.002 0.0006  0.0004   0.06     0.2     1.0  # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
      *      finetune = 1: 5   0.001 0.001   0.001    0.3      0.33    1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
             finetune = 1: 0.00957  0.00712  0.00493  0.00048  0.00746 .01 .01 .01  # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
                           # [9/1/2022] This is what I used in Maxima/Moschata genome paper.

                print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
               burnin = 20000
           * sampfreq = 2
             sampfreq = 10
            * nsample = 100000
              nsample = 40000

            * Threads = 30 1 2
           checkpoint = 10000 10000

*** Note: Make your window wider (140 columns) before running the program.
IMAP
&fileSunhh::write2file($ofn2, $imap_txt, '>>');

&tsmsg("[Msg] The out files are:\n");
&tsmsg("[Msg]   $ofn1 : \n");
&tsmsg("[Msg]   $ofn2 : \n");

