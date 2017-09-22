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
); 
$opts{'cpuN'} //= 10; 

my $help_txt = <<HH; 

perl $0 set02.tab.filt.mjAL.CLV_CLM_CA_CC > set02.tab.filt.mjAL.class_list

-cpuN    [$opts{'cpuN'}]

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


# If 'CC' = 'N', then I may have 'CC_CLM_CLV' and 'CC_CA' together. 
# At least, there is no difference among 'CC CLM CLV' in 'CC_CLM_CLV' type. 
sub class_mjAL {
	my ($clv, $clm, $ca, $cc) = @_; 
	my %h0 = map { $_ => 1 } grep { $_ ne 'N' } ($clv, $clm, $ca, $cc); 
	my @back; 
	for my $b1 (sort keys %h0) {
			
		if ( $b1 eq $cc or $cc eq 'N' ) {
			if ( $b1 eq $ca or $ca eq 'N' ) {
				if ( $b1 eq $clm or $clm eq 'N' ) {
					if ( $b1 eq $clv or $clv eq 'N' ) {
						push(@back, [ $b1, 't5', 'CC_CA_CLM_CLV' ]); 
					} else {
						push(@back, [ $b1, 't3', 'CC_CA_CLM' ]); 
					}
				} elsif ( $b1 eq $clv or $clv eq 'N' ) {
					push(@back, [ $b1, 't4', 'CC_CA_CLV' ]); 
				} else {
					push(@back, [ $b1, 't2', 'CC_CA' ]); 
				}
			} elsif ( $b1 eq $clm or $clm eq 'N' ) {
				if ( $b1 eq $clv or $clv eq 'N' ) {
					push(@back, [ $b1, 't8', 'CC_CLM_CLV' ]); 
				} else {
					push(@back, [ $b1, 't6', 'CC_CLM' ]); 
				}
			} elsif ( $b1 eq $clv or $clv eq 'N' ) {
				push(@back, [ $b1, 't7', 'CC_CLV' ]); 
			} else {
				push(@back, [ $b1, 't1', 'CC' ]); 
			}
		} elsif ( $b1 eq $ca or $ca eq 'N' ) {
			if ( $b1 eq $clm or $clm eq 'N' ) {
				if ( $b1 eq $clv or $clv eq 'N' ) {
					push(@back, [ $b1, 't11', 'CA_CLM_CLV' ]); 
				} else {
					push(@back, [ $b1, 't9', 'CA_CLM' ]); 
				}
			} elsif ( $b1 eq $clv or $clv eq 'N' ) {
				push(@back, [ $b1, 't10', 'CA_CLV' ]); 
			} else {
				push(@back, [ $b1, 't12', 'CA' ]); 
			}
		} elsif ( $b1 eq $clm or $clm eq 'N' ) {
			if ( $b1 eq $clv or $clv eq 'N') {
				push(@back, [ $b1, 't13', 'CLM_CLV' ]); 
			} else {
				push(@back, [ $b1, 't14', 'CLM' ]); 
			}
		} elsif ( $b1 eq $clv ) {
			push(@back, [ $b1, 't15', 'CLV' ]); 
		} else {
			# clv is 'N', this should not happen. 
			&stopErr("[Err] Bad here. $b1 in [$clv, $clm, $ca, $cc]\n"); 
			push(@back, [ $b1, 't0' , 'UN' ]); 
		}
	}
	return(\@back); 
}# class_mjAL() 


