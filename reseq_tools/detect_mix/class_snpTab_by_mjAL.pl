#!/usr/bin/perl -w
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"cpuN:i", 
	"class_mjAL:s", # set02.tab.filt.mjAL.class_list 
	"colN_start:i", # 3
); 
$opts{'cpuN'} //= 10; 
$opts{'colN_start'} //= 3; 


my $help_txt = <<HH; 

perl $0 set02.tab.filt -class_mjAL set02.tab.filt.mjAL.class_list > set02.tab.filt.mjAL_class

 -class_mjAL    [fielname] Required. 
                  Format : 
                  chr             pos     AL      class_1 class_2
                  WM97v2_Chr01    291     T       t5      CC_CA_CLM_CLV
                  WM97v2_Chr01    417     G       t8      CC_CLM_CLV
                  WM97v2_Chr01    417     T       t2      CC_CA

 -colN_start    [$opts{'colN_start'}] 

 -cpuN          [$opts{'cpuN'}]

 -help          [Boolean]
###############################################
 Input example of 'set02.tab.filt.mjAL.class_list' : 
chr             pos     AL      class_1 class_2
WM97v2_Chr01    291     T       t5      CC_CA_CLM_CLV
WM97v2_Chr01    417     G       t8      CC_CLM_CLV
WM97v2_Chr01    417     T       t2      CC_CA
###############################################
 Input example of 'set02.tab.filt' : 
chr             pos     base    indv_base   ind2    ind3    ind4    ind5    ind6
WM97v2_Chr01    291     T       T/T         T/T     T/T     T/T     ./.     T/T
WM97v2_Chr01    417     G       G/G         G/G     G/G     G/G     G/G     G/G
WM97v2_Chr01    432     A       A/A         A/A     A/A     A/A     A/A     A/A
###############################################
 Output format : 
chr             pos     base    sample1         sample2
WM97v2_Chr01    417     T       CC_CA/CC_CA     CC_CLM_CLV/CC_CLM_CLV
###############################################
HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

my @InFp = () ;
if ( !@ARGV ) {
	-t or @InFp = (\*STDIN); 
} else {
	for ( @ARGV ) {
		push( @InFp, &openFH($_,'<') ); 
	}
}

my $pm = &LogInforSunhh::get_pm( $opts{'cpuN'} ); 

my %mjH = &load_mjAL( $opts{'class_mjAL'} ); 

for my $fh ( @InFp ) {
	# header_txt
	my $header_txt = <$fh>; 
	print STDOUT $header_txt; 
	# separate files 
	my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
	my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" ); 
	# process sub-files
	for my $sfn (@sub_fn) {
		my $pid = $pm->start and next;
		open F,'<',"$sfn" or &stopErr("[Err] Failed to open file [$sfn]\n");
		open O,'>',"$sfn.o" or &stopErr("[Err] Failed to open file [$sfn.o]\n");
		while (<F>) {
			chomp; 
			my @ta = &splitL("\t", $_); 
			&class_line( \@ta, \%mjH, $opts{'colN_start'} ); 
			print O join("\t", @ta)."\n"; 
		}
		close O;
		close F;
		$pm->finish;
	}
	$pm->wait_all_children; 
	# merge sub-files
	for my $sfn (@sub_fn) {
		open F,'<',"$sfn.o" or &stopErr("[Err] Failed to open file [$sfn.o]\n");
		while (<F>) {
			print STDOUT $_;
		}
		close F;
	}
	# delete sub-files
	&fileSunhh::_rmtree($wrk_dir); 
}

###############################################
###############################################
sub load_mjAL {
	my $fn = shift; 
	my %back; 
	my $fh = &openFH($fn, '<'); 
	while (<$fh>) {
		chomp; 
		my @ta = &splitL("\t", $_); 
		($ta[0] eq 'chr' or $ta[0] eq '#CHROM') and next; 
		$back{$ta[0]}{$ta[1]}{$ta[2]} = $ta[4]; 
	}
	close($fh); 
	return(%back); 
}# load_mjAL() 

sub class_line {
	my ( $lineAR, $mjHR, $sCN ) = @_; 
	$sCN //= 3; 
	my $chrID = $lineAR->[0]; 
	my $chrP  = $lineAR->[1]; 
	for (my $i=$sCN; $i<@$lineAR; $i++) {
		$lineAR->[$i] =~ m!^([^\s/]+)/([^\s/]+)$! or die "bad genotype [$lineAR->[$i]]\n"; 
		my ($b1, $b2) = ($1, $2); 
		if ( $b1 eq '.' or $b1 eq 'N' ) {
			$lineAR->[$i] = './.'; 
		} else {
			my ($nb1, $nb2) = ('REST', 'REST'); 
			defined $mjHR->{$chrID}{$chrP}{$b1} and $nb1 = $mjHR->{$chrID}{$chrP}{$b1} ; 
			defined $mjHR->{$chrID}{$chrP}{$b2} and $nb2 = $mjHR->{$chrID}{$chrP}{$b2} ; 
			$lineAR->[$i] = "$nb1/$nb2"; 
		}
	}
	return; 
}# class_line () 

