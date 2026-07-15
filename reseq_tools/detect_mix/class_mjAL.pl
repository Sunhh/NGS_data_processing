#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"cpuN:i", 
	"nAsWildcard!",   # legacy: treat a group's 'N' major allele as matching any allele. 
	"onlyDiag!",      # only output diagnostic sites (>= 2 distinct non-N group major alleles). 
); 
$opts{'cpuN'} //= 10; 

my $help_txt = <<HH; 

perl \$0 set02.tab.filt.mjAL.CLV_CLM_CA_CC > set02.tab.filt.mjAL.class_list

Classify each site by which reference groups share the same major allele. Input columns:
  chr  pos  CLV  CLM  CA  CC   (each = that group's major allele, or 'N' if undefined)
Output: chr, pos, AL(allele), class_1(t-code), class_2(group set, e.g. CC_CA / CLM_CLV).

-cpuN          [$opts{'cpuN'}]
-onlyDiag      [Boolean] Skip non-diagnostic sites (fewer than 2 distinct non-N major
               alleles among the groups, i.e. sites that do not distinguish any groups).
-nAsWildcard   [Boolean] LEGACY behaviour: treat a group whose major allele is 'N' as
               matching every allele (so an undefined group is added to the shared class).
               Default (off): 'N' groups are simply left out of the site's class, which
               avoids assigning an undefined group to overlapping classes.
-help

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

# chr           pos     CLV     CLM     CA      CC
# WM97v2_Chr01  291     T       T       N       N
# WM97v2_Chr01  417     G       G       T       N

my @InFp = () ;
if ( !@ARGV ) {
	-t or @InFp = (\*STDIN); 
} else {
	for ( @ARGV ) {
		push( @InFp, &openFH($_,'<') ); 
	}
}

my $pm = &LogInforSunhh::get_pm( $opts{'cpuN'} ); 

my @GRP_ORDER = qw/CC CA CLM CLV/; 
my %CLASS2T = (
	"CC"=>"t1", "CC_CA"=>"t2", "CC_CA_CLM"=>"t3", "CC_CA_CLV"=>"t4", "CC_CA_CLM_CLV"=>"t5", 
	"CC_CLM"=>"t6", "CC_CLV"=>"t7", "CC_CLM_CLV"=>"t8", "CA_CLM"=>"t9", "CA_CLV"=>"t10", 
	"CA_CLM_CLV"=>"t11", "CA"=>"t12", "CLM_CLV"=>"t13", "CLM"=>"t14", "CLV"=>"t15", 
); 

for my $fh ( @InFp ) {
	# header_txt
	my $header_txt = <$fh>; 
	print STDOUT join("\t", qw/chr pos AL class_1 class_2/)."\n"; 
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
			my ($cid, $cP, $clv, $clm, $ca, $cc) = &splitL("\t", $_); 
			my ($cls) = &class_mjAL( $clv, $clm, $ca, $cc ); 
			for my $t1 (@$cls) {
				# @$t1 = ( AL, class_1, class_2 )
				print O join("\t", $cid, $cP, $t1->[0], $t1->[1], $t1->[2])."\n"; 
			}
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

sub class_mjAL {
	my ($clv, $clm, $ca, $cc) = @_; 
	my %al = ( 'CLV'=>$clv, 'CLM'=>$clm, 'CA'=>$ca, 'CC'=>$cc ); 
	my %distinct = map { $_ => 1 } grep { $_ ne 'N' } values %al; 
	my @back; 
	# -onlyDiag : skip sites that do not distinguish any group (< 2 distinct non-N alleles). 
	$opts{'onlyDiag'} and scalar(keys %distinct) < 2 and return(\@back); 
	for my $b1 (sort keys %distinct) {
		my @grps; 
		for my $g (@GRP_ORDER) {
			if ( $al{$g} eq $b1 or ( $opts{'nAsWildcard'} and $al{$g} eq 'N' ) ) {
				push(@grps, $g); 
			}
		}
		my $name  = join('_', @grps); 
		my $tcode = $CLASS2T{$name} // 't0'; 
		push(@back, [ $b1, $tcode, $name ]); 
	}
	return(\@back); 
}# class_mjAL() 
