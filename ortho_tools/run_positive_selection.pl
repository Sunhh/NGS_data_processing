#!/usr/bin/perl
# http://www.ch.embnet.org/CoursEMBnet/PagesPHYL07/Exercises/day2/day2.html
# https://evosite3d.blogspot.com/2011/09/identifying-positive-selection-in.html
### Input tree file '1by1.tree' : ( ( (cucumber.aa.fasta, melon.aa.fasta), watermelon.aa.fasta ), (maxima.aa.fasta, moschata.aa.fasta)#1 );
### Input cds fasta file 'all.cds.fa' : 
### Input protein fasta file 'all.prot.fa' : 
### Input ortholog group file 'OrthologousGroups.csv.CuFam1to2_OGs.1by1' : 
#### OG_ID	cucumber.aa.fasta	melon.aa.fasta	watermelon.aa.fasta	maxima.aa.fasta	moschata.aa.fasta
#### OG0000104	cucumber.Csa4M652860.1	melon.MELO3C017000P1	watermelon.Cla007410	maxima.Cma_004288	moschata.Cmo_006025
#### OG0000177	cucumber.Csa7M045580.1	melon.MELO3C018777P1	watermelon.Cla007714	maxima.Cma_027242	moschata.Cmo_015903
#### OG0000192	cucumber.Csa2M292830.1	melon.MELO3C021359P1	watermelon.Cla010439	maxima.Cma_010246	moschata.Cmo_010352
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 
use Statistics::Distributions; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"in_tree0:s", # 1by1.0.tree
	"in_tree1:s", # 1by1.1.tree
	"in_ctl0:s", # 1by1.0.ctl # H0.ctl : fixed .
	"in_ctl1:s", # 1by1.1.ctl # H1.ctl : alternative 
	"in_cds:s",  # all.cds.fa
	"in_prot:s", # all.prot.fa
	"in_ogs:s",  # OrthologousGroups.csv.CuFam1to2_OGs.1by1
	"repN:i",    # 3 
	"out_test:s", # 
	"exe_muscle:s", # muscle 
	"exe_codeml:s", # codeml
	"keep_tmp!", 
	"cpuN:i",    # 1
	"help!", 
); 
$opts{'exe_muscle'} //= 'muscle'; 
$opts{'exe_codeml'} //= 'codeml'; 
$opts{'cpuN'}       //= 10; 
$opts{'repN'}       //= 3; 

my $help_txt = <<HH;
################################################################################
# perl $0 -in_tree0 1by1.0.tree -in_tree1 1by1.1.tree -in_ctl0 1by1.0.ctl -in_ctl1 1by1.1.ctl -in_cds all.cds.fa -in_prot all.prot.fa -in_ogs og_1by1
#
# -cpuN       [$opts{'cpuN'}]
#
# -out_test   [filename to replace STDOUT]
#
# -repN       [$opts{'repN'}] Times to repeat for each model. Use the maximum lnL among repeats. 
#
# -keep_tmp   [Boolean]
#
# -exe_muscle [$opts{'exe_muscle'}]
# -exe_codeml [$opts{'exe_codeml'}]
#
################################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my $pm = &LogInforSunhh::get_pm($opts{'cpuN'}); 
my $ofh_test = \*STDOUT; 
defined $opts{'out_test'} and $ofh_test = &openFH($opts{'out_test'}, '>'); 

my %cds_h  = %{$fs_obj->save_seq_to_hash( 'faFile'=>$opts{'in_cds'} )}; 
my %prot_h = %{$fs_obj->save_seq_to_hash( 'faFile'=>$opts{'in_prot'} )}; 
defined $opts{'in_tree0'} or die "-in_tree0 required.\n"; 
defined $opts{'in_tree1'} or die "-in_tree1 required.\n"; 
# my %ogs_h  = %{&load_ogs($opts{'in_ogs'})}; 
for (keys %prot_h) { $prot_h{$_}{'seq'} =~ s!\s!!g; }
for (keys %cds_h ) { $cds_h{$_}{'seq'}  =~ s!\s!!g; }

my $ori_dir = &fileSunhh::_abs_path("./"); 
my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
$wrk_dir    = &fileSunhh::_abs_path($wrk_dir); 
&fileSunhh::_copy( $opts{'in_tree0'}, "$wrk_dir/c.0.tree"); 
&fileSunhh::_copy( $opts{'in_tree1'}, "$wrk_dir/c.1.tree"); 
&setup_ctl0("$wrk_dir/c.0.ctl", $opts{'in_ctl0'}); 
&setup_ctl1("$wrk_dir/c.1.ctl", $opts{'in_ctl1'}); 
&tsmsg("[Rec] Go to work dir [$wrk_dir] from [$ori_dir]\n"); 
my @sub_fn = &fileSunhh::dvd_file( $opts{'in_ogs'}, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 1, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" ); 
chdir($wrk_dir); 

for my $sfn (@sub_fn) {
	my $sub_dir = "${sfn}_dir"; 
	$sub_dir = &fileSunhh::_abs_path($sub_dir); 
	mkdir($sub_dir) or &stopErr("[Err] Failed to create sub_dir [$sub_dir]\n"); 
	my $pid = $pm->start and next; 
	my %ogs_h = %{ &load_ogs($sfn) }; 
	&fileSunhh::write2file("$sfn.test", join("\t", qw/OG_ID p_value LRT_delta lnL_1 lnL_0 df Gene_ID tree/)."\n", '>'); 
	&fileSunhh::_copy("$wrk_dir/c.0.ctl", "$sub_dir/"); 
	&fileSunhh::_copy("$wrk_dir/c.1.ctl", "$sub_dir/"); 
	&fileSunhh::_copy("$wrk_dir/c.0.tree",  "$sub_dir/"); 
	&fileSunhh::_copy("$wrk_dir/c.1.tree",  "$sub_dir/"); 
	# goto SUB_END; 
	
	chdir($sub_dir); 
	for my $tr1 ( @{$ogs_h{'OGs'}} ) {
		&tsmsg("[Msg] Processing $tr1->[0] in [$sub_dir]\n"); 
		my $is_err = 0; 
		open OFA_P,'>',"p.fa" or &stopErr("[Err] Failed to write file [$sub_dir/p.fa]\n"); 
		open OFA_C,'>',"c.fa" or &stopErr("[Err] Failed to write file [$sub_dir/c.fa]\n"); 
		for (my $i=1; $i<@$tr1; $i++) {
			print OFA_P ">$ogs_h{'header'}[$i]\n$prot_h{$tr1->[$i]}{'seq'}\n"; 
			print OFA_C ">$ogs_h{'header'}[$i]\n$cds_h{$tr1->[$i]}{'seq'}\n"; 
		}
		close OFA_P; 
		close OFA_C; 
		&exeCmd_1cmd("$opts{'exe_muscle'} -in p.fa -out p.aln -seqtype protein") and do { $is_err = 1; goto OG_END; }; # Default output is FASTA format. 
		my %p_aln = %{ $fs_obj->save_seq_to_hash( 'faFile'=>'p.aln' ) }; 
		my %c_hhh = %{ $fs_obj->save_seq_to_hash( 'faFile'=>'c.fa' ) }; 
		for my $t1 (keys %p_aln) { $p_aln{$t1}{'seq'} =~ s!\s!!g; }
		for my $t1 (keys %c_hhh) { $c_hhh{$t1}{'seq'} =~ s!\s!!g; }
		my $ofh_cALN = &openFH( 'c.aln', '>' ); 
		for my $tk ( sort keys %p_aln ) {
			$p_aln{$tk}{'seq'} =~ s!\s!!g; 
			my $cdsByAA = &fastaSunhh::aa2cds_1seq( $c_hhh{$tk}{'seq'}, $p_aln{$tk}{'seq'} ); 
			unless ( defined $cdsByAA ) {
				$is_err = 1; 
				goto OG_END; 
			}
			print {$ofh_cALN} ">$tk\n$cdsByAA\n"; 
		}
		close($ofh_cALN); 
		my %c_h = %{$fs_obj->save_seq_to_hash( 'faFile'=>'c.aln' )}; 
		my $indv_N = scalar(keys %c_h); 
		my $base_N ; 
		for my $tk (keys %c_h) {
			$c_h{$tk}{'seq'} =~ s!\s!!g; 
			$base_N //= length($c_h{$tk}{'seq'}); 
			$base_N == length($c_h{$tk}{'seq'}) or &stopErr("[Err] Different length for [$tk]\n"); 
		}
		open OPHY,'>',"c.phy" or &stopErr("[Err] Failed to write file [$sub_dir/c.phy]\n"); 
		print OPHY "   $indv_N   $base_N\n"; 
		for my $tk (sort {$c_h{$a}{'Order'} <=> $c_h{$b}{'Order'}} keys %c_h) {
			my $indv_ID = $tk; 
			# $indv_ID    =~ s!\s!_!g; 
			# $indv_ID    =~ tr!)(][:;,!_______!; 
			# length($indv_ID) > 9 and $indv_ID = substr($indv_ID, 0, 9); 
			# $indv_ID = sprintf("%-10s", $indv_ID);
			print OPHY "$indv_ID   $c_h{$tk}{'seq'}\n"; 
		}
		close OPHY; 
		my $ofn_0 = &get_ofn_ctl("./c.0.ctl"); 
		my $ofn_1 = &get_ofn_ctl("./c.1.ctl"); 
		my (%v_0, %v_1); 
		my $cnt_success = 0; 
		for (my $k=0; $k<20; $k++) {
			$cnt_success >= $opts{'repN'} and last; 
			&exeCmd_1cmd("$opts{'exe_codeml'} ./c.0.ctl 1> c.0.out") and next; 
			my %v_0_t   = &infor_mlc($ofn_0); 
			$cnt_success ++; 
			unless ( defined $v_0_t{'lnL'} ) {
				&exeCmd_1cmd("mv $sub_dir ${sub_dir}_bad"); 
				&stopErr("[Err] bad $sfn in [$sub_dir]\n"); 
			}
			&tsmsg("[Msg] Success [$cnt_success] for c.0 of $tr1->[0]\n"); 
			if (defined $v_0{'lnL'}) {
				if ($v_0_t{'lnL'} > $v_0{'lnL'}) {
					%v_0 = %v_0_t; 
				}
			} else {
				%v_0 = %v_0_t; 
			}
		}
		defined $v_0{'lnL'} or do { $is_err = 1; goto OG_END; }; 
		$cnt_success = 0; 
		for (my $k=0; $k<20; $k++) {
			$cnt_success >= $opts{'repN'} and last; 
			&exeCmd_1cmd("$opts{'exe_codeml'} ./c.1.ctl 1> c.1.out") and next; 
			my %v_1_t   = &infor_mlc($ofn_1); 
			$cnt_success ++; 
			&tsmsg("[Msg] Success [$cnt_success] for c.1 of $tr1->[0]\n"); 
			if (defined $v_1{'lnL'}) {
				if ($v_1_t{'lnL'} > $v_1{'lnL'}) {
					%v_1 = %v_1_t; 
				}
			} else {
				%v_1 = %v_1_t; 
			}
		}
		defined $v_1{'lnL'} or do { $is_err = 1; goto OG_END; }; 

		my $delta_LRT = 2 * ($v_1{'lnL'}-$v_0{'lnL'}); 
		my $freedom   = $v_1{'np'}-$v_0{'np'}; 
		&add_ofn(
		  "$sfn.test", 
		  $tr1->[0], 
		  ($delta_LRT < 0) ? 1 : Statistics::Distributions::chisqrprob( $freedom, $delta_LRT ), 
		  $delta_LRT, 
		  $v_1{'lnL'}, 
		  $v_0{'lnL'}, 
		  $freedom, 
		  join(",",@{$tr1}[1 .. $#$tr1]),
		  $v_0{'tree'}
		); 
			
		OG_END: 
		if ( $is_err != 0 ) {
			&add_ofn("$sfn.test", $tr1->[0], ("NA") x 5, join(",",@{$tr1}[1 .. $#$tr1]), 'NA'); 
		}
	}
	chdir($wrk_dir); 

	SUB_END: 
	$pm->finish; 
}
$pm->wait_all_children;
chdir($ori_dir); 

my $has_oh = 0; 
for my $sfn (@sub_fn) {
	open F,'<',"$sfn.test" or &stopErr("[Err] Failed to open file [$sfn.test]\n"); 
	while (<F>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		if ($. == 1) {
			$has_oh == 0 and do { print {$ofh_test} $_."\n";  }; 
			$has_oh = 1; 
			next; 
		}
		print {$ofh_test} "$_\n"; 
	}
	close F; 
}

$opts{'keep_tmp'} or &fileSunhh::_rmtree($wrk_dir); 

sub print_std {
	print STDOUT join("\t", @_)."\n"; 
}
sub add_ofn {
	my $ofn = shift; 
	&fileSunhh::write2file($ofn,join("\t", @_)."\n",'>>'); 
}

sub infor_mlc {
	my $fn = shift; 
	my %back; 
	my $fh = &openFH($fn, '<'); 
	while (<$fh>) {
		if (m/^lnL\(ntime:\s*(\d+)\s*np:\s*(\d+)\):\s*(\S+)/) {
			$back{'ntime'} = $1; 
			$back{'np'}    = $2; 
			$back{'lnL'}   = $3; 
		} elsif (m/^Bayes Empirical Bayes/) {
			<$fh>; 
			while (my $tl = <$fh>) {
				$tl =~ m/^$/ and last; 
				$tl =~ m/^\s+/ or last; 
				$tl =~ m/^\s+(\d+)\s+(\S+)\s+(\S+)$/ or &stopErr("[Err] Failed to match BEB on in file [$fn] : $tl\n"); 
				push(@{$back{'BEB'}}, join(":", $1,$2,$3)); 
			}
		} elsif (m/^tree\s+length\s+=/) {
			<$fh>; <$fh>; <$fh>; 
			my $tree = <$fh>; 
			chomp($tree); 
			$back{'tree'} = $tree; 
		}
	}
	close($fh); 
	return(%back); 
}# infor_mlc() 

sub get_ofn_ctl {
	my $fn = shift; 
	my $fh = &openFH($fn,'<'); 
	while (<$fh>) {
		m/^\s*outfile\s*=\s*(\S+)/ or next; 
		return $1; 
		last; 
	}
	close ($fh); 
}# get_ofn_ctl () 

sub load_ogs {
	my $fn = shift; 
	my $fh = &openFH($fn,'<'); 
	my %back; 
	my $has_header = 0; 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		if ( $has_header == 0 ) {
			$has_header = 1; 
			$back{'header'} = [ @ta ]; 
			next; 
		}
		push(@{$back{'OGs'}}, [@ta]); 
	}
	close($fh); 
	return(\%back); 
}# load_ogs () 


sub setup_ctl0 {
	my ($ofn, $ifn) = @_; 
	my $ofh = &openFH($ofn, '>'); 
	if ( defined $ifn ) {

my $ifh = &openFH($ifn, '<'); 
while (<$ifh>) {
	chomp; 
	s!^(\s*seqfile\s*=\s*)(\S+)(\s*)(\*|$)!$1c.phy$3$4!; 
	s!^(\s*treefile\s*=\s*)(\S+)(\s*)(\*|$)!$1c.0.tree$3$4!; 
	s!^(\s*outfile\s*=\s*)(\S+)(\s*)(\*|$)!$1c.0.mlc$3$4!; 
	print {$ofh} "$_\n"; 
}
close ($ifh); 

	} else {
# This is for branch-site model referring: 
#    'http://www.ch.embnet.org/CoursEMBnet/PagesPHYL07/Exercises/day2/lysozyme_M0.ctl'
#    'examples/lysozyme/lysozymeLarge.ctl'
# This is slightly different from 'https://evosite3d.blogspot.com/2011/09/identifying-positive-selection-in.html'
print {$ofh} <<CTL0; 
seqfile  = c.phy              * sequence data file name
treefile = c.0.tree       * tree structure file name
outfile  = c.0.mlc  * main result file name

  noisy = 3     * 0,1,2,3,9: how much rubbish on the screen
verbose = 0     * 1: detailed output, 0: concise output
runmode = 0     * 0: user tree;  1: semi-automatic;  2: automatic
                * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

  seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
    clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
   aaDist = 0   * 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
    model = 2   * models for codons:
                * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
  NSsites = 2   * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete;
                * 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal
    icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
*    Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
    kappa = 3   * initial or fixed kappa
fix_omega = 1   * 1: omega or omega_1 fixed, 0: estimate
    omega = 1   * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0  * different alphas for genes
        ncatG = 10  * # of categories in dG of NSsites models


       getSE = 0       * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0       * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
  Small_Diff = .5e-6   * Default value.
*   cleandata = 1       * remove sites with ambiguity data (1:yes, 0:no)?
*       method = 1   * 0: simultaneous; 1: one branch at a time
* fix_blength = 0       * 0: ignore, -1: random, 1: initial, 2: fixed 
CTL0

	}
	close ($ofh); 
}# setup_ctl0() 

sub setup_ctl1 {
	my ($ofn, $ifn) = @_; 
	my $ofh = &openFH($ofn, '>'); 
	if ( defined $ifn ) {

my $ifh = &openFH($ifn, '<'); 
while (<$ifh>) {
	chomp; 
	s!^(\s*seqfile\s*=\s*)(\S+)(\s*)(\*|$)!$1c.phy$3$4!; 
	s!^(\s*treefile\s*=\s*)(\S+)(\s*)(\*|$)!$1c.1.tree$3$4!; 
	s!^(\s*outfile\s*=\s*)(\S+)(\s*)(\*|$)!$1c.1.mlc$3$4!; 
	print {$ofh} "$_\n"; 
}
close ($ifh); 

	} else {

print {$ofh} <<CTL1; 
seqfile  = c.phy              * sequence data file name
treefile = c.1.tree       * tree structure file name
outfile  = c.1.mlc  * main result file name

  noisy = 3     * 0,1,2,3,9: how much rubbish on the screen
verbose = 0     * 1: detailed output, 0: concise output
runmode = 0     * 0: user tree;  1: semi-automatic;  2: automatic
                * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

  seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
    clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
   aaDist = 0   * 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
    model = 2   * models for codons:
                * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
  NSsites = 2   * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete;
                * 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal
    icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
*    Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
    kappa = 3   * initial or fixed kappa
fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate
    omega = 1   * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0  * different alphas for genes
        ncatG = 10  * # of categories in dG of NSsites models


       getSE = 0       * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0       * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
  Small_Diff = .5e-6   * Default value.
*   cleandata = 1       * remove sites with ambiguity data (1:yes, 0:no)?
*       method = 1   * 0: simultaneous; 1: one branch at a time
* fix_blength = 0       * 0: ignore, -1: random, 1: initial, 2: fixed
CTL1

	}
	close ($ofh); 
}# setup_ctl1() 

