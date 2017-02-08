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
	"in_tree:s", # 1by1.tree
	"in_ctl0:s", # 1by1.0.ctl # fixed
	"in_ctl1:s", # 1by1.1.ctl # alternative 
	"in_cds:s",  # all.cds.fa
	"in_prot:s", # all.prot.fa
	"in_ogs:s",  # OrthologousGroups.csv.CuFam1to2_OGs.1by1
	"exe_muscle:s", # muscle 
	"exe_codeml:s", # codeml
	"help!", 
); 
$opts{'exe_muscle'} //= 'muscle'; 
$opts{'exe_codeml'} //= 'codeml'; 

my $help_txt = <<HH;
################################################################################
# perl $0 -in_tree 1by1.tree -in_ctl0 1by1.0.ctl -in_ctl1 1by1.1.ctl -in_cds all.cds.fa -in_prot all.prot.fa -in_ogs og_1by1
#
#
################################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my %cds_h  = %{$fs_obj->save_seq_to_hash( 'faFile'=>$opts{'in_cds'} )}; 
my %prot_h = %{$fs_obj->save_seq_to_hash( 'faFile'=>$opts{'in_prot'} )}; 
defined $opts{'in_tree'} or die "-in_tree required.\n"; 
my %ogs_h  = %{&load_ogs($opts{'in_ogs'})}; 
for (keys %prot_h) { $prot_h{$_}{'seq'} =~ s!\s!!g; }
for (keys %cds_h ) { $cds_h{$_}{'seq'}  =~ s!\s!!g; }

my $ori_dir = &fileSunhh::_abs_path("./"); 
my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
&fileSunhh::_copy( $opts{'in_tree'}, "$wrk_dir/c.tree"); 
&setup_ctl0("$wrk_dir/c.0.ctl", $opts{'in_ctl0'}); 
&setup_ctl1("$wrk_dir/c.1.ctl", $opts{'in_ctl1'}); 


print STDOUT join("\t", qw/OG_ID p_value LRT_delta lnL_1 lnL_0 df Gene_ID tree/)."\n"; 
&tsmsg("[Rec] Go to work dir [$wrk_dir] from [$ori_dir]\n"); 
chdir($wrk_dir); 
OG: 
for my $tr1 (@{$ogs_h{'OGs'}}) {
	open OFA_P,'>',"p.fa" or die; 
	open OFA_C,'>',"c.fa" or die; 
	for (my $i=1; $i<@$tr1; $i++) {
		print OFA_P ">$ogs_h{'header'}[$i]\n$prot_h{$tr1->[$i]}{'seq'}\n"; 
		print OFA_C ">$ogs_h{'header'}[$i]\n$cds_h{$tr1->[$i]}{'seq'}\n"; 
	}
	close OFA_C; 
	close OFA_P; 
	&exeCmd_1cmd("$opts{'exe_muscle'} -in p.fa -out p.aln -seqtype protein") and do { &print_std($tr1->[0], ("NA") x 5, join(",",@{$tr1}[1 .. $#$tr1]), 'NA'); next; }; # Default output is FASTA format. 
	my %p_aln = %{ $fs_obj->save_seq_to_hash( 'faFile'=>'p.aln' ) }; 
	my %c_hhh = %{ $fs_obj->save_seq_to_hash( 'faFile'=>'c.fa' ) }; 
	for my $t1 (keys %p_aln) { $p_aln{$t1}{'seq'} =~ s!\s!!g; }
	for my $t1 (keys %c_hhh) { $c_hhh{$t1}{'seq'} =~ s!\s!!g; }
	my $ofh_cALN = &openFH( 'c.aln', '>' ); 
	for my $tk ( sort keys %p_aln ) {
		$p_aln{$tk}{'seq'} =~ s!\s!!g; 
		my $cdsByAA = &fastaSunhh::aa2cds_1seq( $c_hhh{$tk}{'seq'}, $p_aln{$tk}{'seq'} ); 
		unless ( defined $cdsByAA ) {
			&print_std($tr1->[0], ("NA") x 5, join(",",@{$tr1}[1 .. $#$tr1]), 'NA'); 
			next OG; 
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
		$base_N == length($c_h{$tk}{'seq'}) or die "Different length for [$tk]\n"; 
	}
	open OPHY,'>',"c.phy" or die; 
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
	&exeCmd_1cmd("$opts{'exe_codeml'} ./c.0.ctl 1> c.0.out") and do { &print_std($tr1->[0], ("NA") x 5, join(",",@{$tr1}[1 .. $#$tr1]), 'NA'); next; }; 
	&exeCmd_1cmd("$opts{'exe_codeml'} ./c.1.ctl 1> c.1.out") and do { &print_std($tr1->[0], ("NA") x 5, join(",",@{$tr1}[1 .. $#$tr1]), 'NA'); next; }; 
	my $ofn_0 = &get_ofn_ctl("./c.0.ctl"); 
	my $ofn_1 = &get_ofn_ctl("./c.1.ctl"); 
	my %v_0   = &infor_mlc($ofn_0); 
	my %v_1   = &infor_mlc($ofn_1); 
	# &exeCmd_1cmd("rm -f p.fa p.aln c.fa c.aln c.phy $ofn_0 $ofn_1"); 
	my $delta_LRT = 2 * ($v_1{'lnL'}-$v_0{'lnL'}); 
	my $freedom   = $v_1{'np'}-$v_0{'np'}; 
	&print_std(
	  $tr1->[0], 
	  Statistics::Distributions::chisqrprob( $freedom, $delta_LRT ), 
	  $delta_LRT, 
	  $v_1{'lnL'}, 
	  $v_0{'lnL'}, 
	  $freedom, 
	  join(",",@{$tr1}[1 .. $#$tr1]),
	  $v_0{'tree'}
	); 
}
chdir($ori_dir); 
&fileSunhh::_rmtree($wrk_dir); 

sub print_std {
	print STDOUT join("\t", @_)."\n"; 
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
				$tl =~ m/^\s+(\d+)\s+(\S+)\s+(\S+)$/ or die "$tl\n"; 
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
	s!^(\s*treefile\s*=\s*)(\S+)(\s*)(\*|$)!$1c.tree$3$4!; 
	s!^(\s*outfile\s*=\s*)(\S+)(\s*)(\*|$)!$1c.0.mlc$3$4!; 
	print {$ofh} "$_\n"; 
}
close ($ifh); 

	} else {

print {$ofh} <<CTL0; 
seqfile  = c.phy              * sequence data file name
treefile = c.tree       * tree structure file name
outfile  = c.0.mlc  * main result file name

  noisy = 9     * 0,1,2,3,9: how much rubbish on the screen
verbose = 1     * 1: detailed output, 0: concise output
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
    Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
    kappa = 2   * initial or fixed kappa
fix_omega = 1   * 1: omega or omega_1 fixed, 0: estimate
    omega = 1   * initial or fixed omega, for codons or codon-based AAs

       getSE = 0       * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0       * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
  Small_Diff = .45e-6  * Default value.
   cleandata = 1       * remove sites with ambiguity data (1:yes, 0:no)?
 fix_blength = 0       * 0: ignore, -1: random, 1: initial, 2: fixed 
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
	s!^(\s*treefile\s*=\s*)(\S+)(\s*)(\*|$)!$1c.tree$3$4!; 
	s!^(\s*outfile\s*=\s*)(\S+)(\s*)(\*|$)!$1c.1.mlc$3$4!; 
	print {$ofh} "$_\n"; 
}
close ($ifh); 

	} else {

print {$ofh} <<CTL1; 
seqfile  = c.phy              * sequence data file name
treefile = c.tree       * tree structure file name
outfile  = c.1.mlc        * main result file name

  noisy = 9     * 0,1,2,3,9: how much rubbish on the screen
verbose = 1     * 1: detailed output, 0: concise output
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
    Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
    kappa = 2   * initial or fixed kappa
fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate
    omega = 1   * initial or fixed omega, for codons or codon-based AAs

       getSE = 0       * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0       * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
  Small_Diff = .45e-6  * Default value.
   cleandata = 1       * remove sites with ambiguity data (1:yes, 0:no)?
 fix_blength = 0       * 0: ignore, -1: random, 1: initial, 2: fixed 
CTL1

	}
	close ($ofh); 
}# setup_ctl1() 

