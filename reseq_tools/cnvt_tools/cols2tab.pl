#!/usr/bin/perl 
# 2016-12-01 Add -cpuN 
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use SNP_tbl; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"cpuN:i", 
	"colN_start:i", # 3
	"colN_ref:i",   # '' 
	"noSrtAllele!", # Do not sort alleles within genotype by character order. 
	"noHeader!",    # 
); 

$opts{'cpuN'} //= 1; 
$opts{'cpuN'} = int($opts{'cpuN'}); 
$opts{'cpuN'} < 1 and $opts{'cpuN'} = 1; 
$opts{'colN_start'} //= 3; 


my $help_txt = <<HH; 

perl $0 in_sMao.snp.cols > in_sMao.snp.tab

-help
-cpuN           [$opts{'cpuN'}] 

-colN_start     [$opts{'colN_start'}]
-colN_ref       [] If this is used, the intact column will be copied just before colN_start. 

-noHeader       [Boolean] If given, the 1st line will be modified too. 

-noSrtAllele    [Boolean] If given, I won't sort alleles within genotype, 
                          so the 'AT' will be 'A/T', and 'TA' will be 'T/A'. 
                          If not given, 'AT' or 'TA' will always be 'A/T'. 

Format of in_sMao.snp.cols : 
  chr \\t pos \\t base(ref) \\t sample1 \\t sample2 ... 
  c1  \\t 100 \\t A         \\t AT      \\t G
  c2  \\t 10  \\t G         \\t N       \\t *
  c3  \\t 5   \\t T         \\t A+ATG   \\t *A+ATG
  c4  \\t 8   \\t N         \\t A*      \\t R
  c4  \\t 10  \\t A         \\t AGT     \\t H

Format of vcf.tab : 
  chr \\t pos \\t base(ref) \\t sample1    \\t sample2 ...
  c1  \\t 100 \\t A         \\t A/T        \\t G/G
  c2  \\t 10  \\t G         \\t ./.        \\t */*
  c3  \\t 5   \\t T         \\t AATG/AATG  \\t AATG/*
  c4  \\t 8   \\t N         \\t A/*        \\t A/G
  c4  \\t 10  \\t A         \\t ./.        \\t ./.

By default, I'll keep column 'base' as single letter, and 

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

my @InFp = (); 
if ( !@ARGV ) {
	-t or @InFp = (\*STDIN); 
} else {
	for (@ARGV) {
		push( @InFp, &openFH($_, '<') ); 
	}
}
my $pm; 
$opts{'cpuN'} > 1 and $pm = &LogInforSunhh::get_pm( $opts{'cpuN'} ); 

my %d2b_list = &SNP_tbl::get_diploid_d2b(); 

my %used; 
$used{'bad_geno'} = {}; 

for my $fh ( @InFp ) {
	if ( defined $pm ) {
		# In this case, %{$used{'bad_geno'}} is not completely recorded. 
		$opts{'cpuN'} > 1 or &stopErr("[Err] cpuN not > 1\n"); 
		# Keep header 
		my $header_txt = ''; 
		unless ( $opts{'noHeader'} ) {
			$header_txt = <$fh>; 
			chomp($header_txt); 
			my @ta = &splitL("\t", $header_txt); 
			my @med; 
			defined $opts{'colN_ref'} and @med = ($ta[$opts{'colN_ref'}]); 
			my @tb = @ta[ $opts{'colN_start'} .. $#ta ]; 
			$header_txt = join("\t", @ta[ 0 .. ($opts{'colN_start'} - 1) ], @med, @tb)."\n";  
		}
		# separate files 
		my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
		my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" ); 
		# process sub-files
		for my $sfn (@sub_fn) {
			my $pid = $pm->start and next;
			open F,'<',"$sfn" or die; 
			open O,'>',"$sfn.o" or die; 
			while (<F>) {
				chomp; 
				my @ta = &splitL( "\t", $_ ); 
				my @med; 
				defined $opts{'colN_ref'} and @med = ($ta[$opts{'colN_ref'}]); 
				my @tb = @ta[ $opts{'colN_start'} .. $#ta ]; 
				&aref_cols2tab( \@tb, $used{'bad_geno'}, \%d2b_list, $opts{'noSrtAllele'} ); 
				print O join("\t", @ta[ 0 .. ($opts{'colN_start'} - 1) ], @med, @tb)."\n";  
			}
			close O; 
			close F; 
			$pm->finish; 
		}
		$pm->wait_all_children; 
		# merge sub-results
		print STDOUT $header_txt; 
		for my $sfn (@sub_fn) {
			open F,'<',"$sfn.o" or die; 
			while ( <F> ) {
				print STDOUT $_; 
			}
			close F; 
		}
		# delete sub-files
		&fileSunhh::_rmtree($wrk_dir); 
	} else {
		my $curr_lineN = 0; 
		while (<$fh>) { 
			$curr_lineN ++; 
			chomp; 
			my @ta=split(/\t/, $_); 
			my @med; 
			defined $opts{'colN_ref'} and @med = ($ta[$opts{'colN_ref'}]); 
			my @tb = @ta[ $opts{'colN_start'} .. $#ta ]; 

			if ( $curr_lineN == 1 and !$opts{'noHeader'} ) {
				print STDOUT join("\t", @ta[ 0 .. ($opts{'colN_start'} - 1) ], @med, @tb)."\n"; 
				next; 
			}
			&aref_cols2tab( \@tb, $used{'bad_geno'}, \%d2b_list, $opts{'noSrtAllele'} ); 
			print STDOUT join("\t", @ta[ 0 .. ($opts{'colN_start'} - 1) ], @med, @tb)."\n";  
		}
	}
	close($fh); 
}

=head1 aref_cols2tab( \@input_cols, \%used_bad_genotype, $doNotSortAlleles, \%d2b_list )

 $doNotSortAlleles is 0 by default, which means sort alleles (FALSE). 

 @input_cols and %used_bad_genotype will be modified. 

 %d2b_list comes from &SNP_tbl::get_diploid_d2b(); 

 Input array 
  @input_cols = ( 'n',   'a',   'A',   '*',   '-',   'W',   'TA',  '*T', 'T+AA',    'TT+AA', '*T+AA', 'GAT', 'B'  ) 
 will be changed to 
  @input_cols = ( './.', 'A/A', 'A/A', '*/*', '*/*', 'A/T', 'A/T', 'T/*','TAA/TAA', 'T/TAA', '*/TAA', './.', './.')
 Input %used_bad_genotype will have 
  $used_bad_genotype{'GAT'} = 1; 
  $used_bad_genotype{'B'} = 1; 

 If $doNotSortAlleles is 1, 'TA' will be 'T/A' instead of 'A/T'; 

Return : undef(); 


=cut
sub aref_cols2tab {
	my ($ar, $href_used, $href_d2b, $noSrtAl) = @_; 
	defined $href_used or &stopErr("[Err] Need \%bad_geno_record\n"); 
	defined $href_d2b or &stopErr("[Err] Need \%d2b_list\n"); 
	$noSrtAl //= 0; 
	for my $tb (@$ar) {
		$tb = uc($tb); 
		$tb eq '-' and $tb = '*'; 
		my @tc; 
		if ( $tb =~ m!^([ATGC\*])$! ) {
			@tc = ($1, $1); 
		} elsif ( $tb eq 'N' ) {
			@tc = ('.', '.'); 
		} elsif ( defined $href_d2b->{$tb} ) {
			@tc = @{$href_d2b->{$tb}}; 
			# if ( $noSrtAl ) {
			# } else {
			#	@tc = sort @{$href_d2b->{$tb}}; 
			#	$tc[0] eq '*' and @tc[0,1] = @tc[1,0]; 
			#}
		} elsif ( $tb =~ m!^([ATGC\*])([ATGC\*])$! ) {
			@tc = ($1,$2); 
			unless ( $noSrtAl ) {
				@tc = sort @tc; 
				$tc[0] eq '*' and @tc[0,1] = @tc[1,0]; 
			}
		} elsif ( $tb =~ m!^([ATGC])\+([ATGCN]+)$! ) {
			@tc = ("$1$2", "$1$2"); 
		} elsif ( $tb =~ m!^([ATGC\*])([ATGC])\+([ATGC]+)$! ) {
			@tc = ($1, "$2$3"); 
		} elsif ( $tb =~ m!^\+! ) {
			@tc = ('.', '.'); 
		} else {
			unless (defined $href_used->{$tb}) {
				$href_used->{$tb} = 1; 
				&tsmsg("[Wrn] Ignore bad genotype [$tb]\n"); 
			}
			@tc = ('.', '.'); 
		}
		$tb = $tc[0] . '/' . $tc[1]; 
	}
	return; 
}# aref_cols2tab() 

