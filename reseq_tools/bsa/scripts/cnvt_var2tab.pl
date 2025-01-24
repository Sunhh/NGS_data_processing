#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"highID:s@", "lowID:s@", 
	"highParentID:s", 
	"lowParentID:s", 
	"minRefAF:f", "maxRefAF:f", # Min/Max reference allele frequency accepted; 
	"minTotalDepth:i", "maxTotalDepth:i", # Min/Max total depth accepted; 
	"minSampleDepth:i", # accepted: Minimum depth of each sample compared 
	"minGQ:f", # accepted: min Genotype quality of each sample 
	"onlyBiAllele!", 
); 

$opts{'minRefAF'}       //= -1; 
$opts{'maxRefAF'}       //= -1; 
$opts{'minTotalDepth'}  //= -1; 
$opts{'maxTotalDepth'}  //= -1; 
$opts{'minSampleDepth'} //= -1; 
$opts{'minGQ'}          //= -1; 

my $help_txt = <<HH; 
################################################################################
# perl $0 input.VariantsToTable.table -highID F2BMut -lowID F2H > input.bsaTab
#
# -help 
#
# -highParentID    [sampleID] homozygous high parent; I always assume highParent is ALT; 
# -lowParentID     [sampleID] homozygous low parent;  I always assume lowParent  is REF; 
#
# -minRefAF        [0-1]
# -maxRefAF        [0-1]
# -minTotalDepth   [Num] Any sites with total depth ==0 will be removed. 
# -maxTotalDepth   [Num]
# -minSampleDepth  [Num]
# -minGQ           [float]
#
# -onlyBiAllele    [Boolean] Ignore sites with more than one allele in ALT. 
################################################################################
# Data with NA will be skipped. 
# Output file format : 
#   CHROM   POS   n2(Ref_AD.high) n4(ALT_AD.high) n1(REF_AD.low) n3(ALT_AD.low)
#   WK1_R1  73    7               6               16             8
#   ....
################################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
defined $opts{'highID'} or &LogInforSunhh::usage($help_txt); 
defined $opts{'lowID'} or &LogInforSunhh::usage($help_txt); 

my %gg; 
&set_Glob(); 
sub set_Glob {
	for (@{$opts{'highID'}}) { $gg{'highID'}{$_} = 1; }
	for (@{$opts{'lowID'}})  { $gg{'lowID'}{$_} = 1; }
	if ( defined $opts{'highParentID'} ) {
		$gg{'highParentID'}{$opts{'highParentID'}} = 1; 
	}
	if ( defined $opts{'lowParentID'} ) {
		$gg{'lowParentID'}{$opts{'lowParentID'}} = 1; 
	}
}# set_Glob() 

my %colN; 
SITE: 
while (<>) {
	chomp; 
	my @ta = &splitL("\t", $_); 
	if ( $. == 1 ) {
		for (my $i=0; $i<@ta; $i++) {
			if      ($ta[$i] =~ m!^(\S+)\.AD$!) {
				if      ( defined $gg{'highID'}{$1} ) {
					push(@{$colN{'high_AD'}}, $i); 
				} elsif ( defined $gg{'lowID'}{$1} ) {
					push(@{$colN{'low_AD'}}, $i); 
				} elsif ( defined $gg{'highParentID'}{$1} ) {
					push(@{$colN{'highParent_AD'}}, $i); 
				} elsif ( defined $gg{'lowParentID'}{$1} ) {
					push(@{$colN{'lowParent_AD'}}, $i); 
				}
			} elsif ($ta[$i] =~ m!^(\S+)\.DP$!) {
				if      ( defined $gg{'highID'}{$1} ) {
					push(@{$colN{'high_DP'}}, $i); 
				} elsif ( defined $gg{'lowID'}{$1} ) {
					push(@{$colN{'low_DP'}}, $i); 
				} elsif ( defined $gg{'highParentID'}{$1} ) {
					push(@{$colN{'highParent_DP'}}, $i); 
				} elsif ( defined $gg{'lowParentID'}{$1} ) {
					push(@{$colN{'lowParent_DP'}}, $i); 
				}
			} elsif ($ta[$i] =~ m!^(\S+)\.GQ$!) {
				if      ( defined $gg{'highID'}{$1} ) {
					push(@{$colN{'high_GQ'}}, $i); 
				} elsif ( defined $gg{'lowID'}{$1} ) {
					push(@{$colN{'low_GQ'}}, $i); 
				} elsif ( defined $gg{'highParentID'}{$1} ) {
					push(@{$colN{'highParent_GQ'}}, $i); 
				} elsif ( defined $gg{'lowParentID'}{$1} ) {
					push(@{$colN{'lowParent_GQ'}}, $i); 
				}
			} elsif ($ta[$i] =~ m!^(\S+)\.PL$!) {
				if      ( defined $gg{'highID'}{$1} ) {
					push(@{$colN{'high_PL'}}, $i); 
				} elsif ( defined $gg{'lowID'}{$1} ) {
					push(@{$colN{'low_PL'}}, $i); 
				} elsif ( defined $gg{'highParentID'}{$1} ) {
					push(@{$colN{'highParent_PL'}}, $i); 
				} elsif ( defined $gg{'lowParentID'}{$1} ) {
					push(@{$colN{'lowParent_PL'}}, $i); 
				}
			} elsif ($ta[$i] =~ m!^CHROM$!) {
				$colN{'CHROM'} = [$i]; 
			} elsif ($ta[$i] =~ m!^POS$!) {
				$colN{'POS'} = [$i]; 
			} else {
				&tsmsg("[Wrn] Skip col-$i of [$ta[$i]]\n"); 
			}
		}
		for my $k1 (qw/high_AD high_DP low_AD low_DP CHROM POS/) {
			( defined $colN{$k1} and @{$colN{$k1}} > 0 ) or &stopErr("[Err] failed to find column for [$k1]\n"); 
		}
		if (defined $opts{'highParentID'}) {
			defined $colN{'highParent_AD'} or &stopErr("[Err] Failed to find column for highParent_AD\n"); 
		}
		if (defined $opts{'lowParentID'}) {
			defined $colN{'lowParent_AD'} or &stopErr("[Err] Failed to find column for lowParent_AD\n"); 
		}
		print STDOUT join("\t", qw/CHROM POS n2_RefAD_H n4_AltAD_H n1_RefAD_L n3_AltAD_L/)."\n"; 
		next SITE; 
	}
	if ($opts{'onlyBiAllele'}) {
		$ta[ $colN{'high_AD'}[0] ] =~ m!^\d+,\d+$! or next SITE; 
		$ta[ $colN{'low_AD'}[0]  ] =~ m!^\d+,\d+$! or next SITE; 
	}
	# Set each value for REF/ALT; 
	my %curr; 
	### Set alt-high and ref-low column in AD list; 
	my ($altHighC, $refLowC) = (1,0); 
	if ( defined $opts{'lowParentID'} ) {
		my @t1 = split(/,/, $ta[ $colN{'lowParent_AD'}[0] ]); 
		my ($max_i, $max_v); 
		for (my $i1 = 0; $i1<@t1; $i1++) {
			$t1[$i1] > 0 or next; 
			$max_i //= $i1; 
			$max_v //= $t1[$i1]; 
			if ($max_v < $t1[$i1]) {
				$max_i = $i1; 
				$max_v = $t1[$i1]; 
			}
		}
		defined $max_i or do { &tsmsg("[Wrn] Skip line with bad lowParent_AD [$ta[ $colN{'lowParent_AD'}[0] ]] at [@ta[0,1]]\n"); next SITE; }; 
		$refLowC = $max_i; 
	}
	if ( defined $opts{'highParentID'} ) {
		my @t2 = split(/,/, $ta[ $colN{'highParent_AD'}[0] ]); 
		my ($max_i2, $max_v2); 
		for (my $i2=0; $i2<@t2; $i2++) {
			$t2[$i2] > 0 or next; 
			$max_i2 //= $i2; 
			$max_v2 //= $t2[$i2]; 
			if ($max_v2 < $t2[$i2]) {
				$max_i2 = $i2; 
				$max_v2 = $t2[$i2]; 
			}
		}
		defined $max_i2 or do { &tsmsg("[Wrn] Skip line with bad highParent_AD [$ta[ $colN{'lowParent_AD'}[0] ]] at [@ta[0,1]]\n"); next SITE; }; 
		$altHighC = $max_i2; 
	}
	if ( defined $opts{'lowParentID'} ) {
		if ( defined $opts{'highParentID'} ) {
			if ($refLowC == $altHighC) {
				&tsmsg("[Err] lowP_AD =[$ta[ $colN{'lowParent_AD'}[0] ]];\n"); 
				&tsmsg("[Err] highP_AD=[$ta[ $colN{'highParent_AD'}[0] ]];\n"); 
				&tsmsg("[Err] Skip bad line for high and low parent IDs genotype: $_\n"); 
				next SITE; 
			}
		} else {
			my @t1 = split(/,/, $ta[ $colN{'lowParent_AD'}[0] ]); 
			for (my $i1=0; $i1<@t1; $i1++) {
				$i1 == $refLowC and next; 
				$altHighC = $i1; 
				last; 
			}
		}
	} elsif ( defined $opts{'highParentID'} ) {
		my @t2 = split(/,/, $ta[ $colN{'highParent_AD'}[0] ]); 
		for (my $i2=0; $i2<@t2; $i2++) {
			$i2 == $altHighC and next; 
			$refLowC = $i2; 
			last; 
		}
	} else {
		; 
	}
	### count AD for ref-low and alt-high; 
	for my $c1 (@{$colN{'low_AD'}}) {
		my @t1 = split(/,/, $ta[$c1]); 
		$t1[ $refLowC  ] =~ s!^\s+|\s+$!!g; $t1[ $refLowC ]  =~ s!^NA$!0!i; 
		$t1[ $altHighC ] =~ s!^\s+|\s+$!!g; $t1[ $altHighC ] =~ s!^NA$!0!i; 
		$curr{'ref_AD_L'} += $t1[$refLowC]; # previous ref_AD_H; 
		$curr{'alt_AD_L'} += $t1[$altHighC]; 
	}
	for my $c1 (@{$colN{'high_AD'}}) {
		my @t1 = split(/,/, $ta[$c1]); 
		$t1[$refLowC]  =~ s!^\s+|\s+$!!g; $t1[$refLowC]  =~ s!^NA$!0!i; 
		$t1[$altHighC] =~ s!^\s+|\s+$!!g; $t1[$altHighC] =~ s!^NA$!0!i; 
		$curr{'ref_AD_H'} += $t1[$refLowC]; 
		$curr{'alt_AD_H'} += $t1[$altHighC]; 
	}
	$curr{'totalDep'} = $curr{'ref_AD_H'}+$curr{'alt_AD_H'}+$curr{'ref_AD_L'}+$curr{'alt_AD_L'}; 
	$curr{'sampleDep_H'} = $curr{'ref_AD_H'}+$curr{'alt_AD_H'}; 
	$curr{'sampleDep_L'} = $curr{'ref_AD_L'}+$curr{'alt_AD_L'}; 
	$curr{'sampleDep_H'} > 0 or next SITE; 
	$curr{'sampleDep_L'} > 0 or next SITE; 
	if ($opts{'minSampleDepth'} > 0) {
		$curr{'sampleDep_H'} >= $opts{'minSampleDepth'} or next SITE; 
		$curr{'sampleDep_L'} >= $opts{'minSampleDepth'} or next SITE; 
	}
	$curr{'totalDep'} > 0 or next SITE; 
	$opts{'minTotalDepth'} > 0 and $curr{'totalDep'} < $opts{'minTotalDepth'} and next SITE; 
	$opts{'maxTotalDepth'} > 0 and $curr{'totalDep'} > $opts{'maxTotalDepth'} and next SITE; 
	$curr{'refAF'} = ($curr{'ref_AD_H'}+$curr{'ref_AD_L'}) / $curr{'totalDep'}; # Allele frequency of reference allele in both bulks; 
	$opts{'minRefAF'} > 0 and $curr{'refAF'} < $opts{'minRefAF'} and next SITE; 
	$opts{'maxRefAF'} > 0 and $curr{'refAF'} > $opts{'maxRefAF'} and next SITE; 

	if ($opts{'minGQ'} > -1) {
		if ( defined $colN{'high_GQ'} and @{$colN{'high_GQ'}} == 1 ) {
			$ta[ $colN{'high_GQ'}[0] ] eq 'NA' and $ta[ $colN{'high_GQ'}[0] ] = 0; 
			$ta[ $colN{'high_GQ'}[0] ] >= $opts{'minGQ'} or next SITE; 
		} elsif ( defined $colN{'high_GQ'} and @{$colN{'high_GQ'}} != 1 ) {
			; 
		} else {
			&tsmsg("[Wrn] Failed to find high_GQ for site @ta[ $colN{'CHROM'}[0], $colN{'POS'}[0] ]\n"); 
		}
		if ( defined $colN{'low_GQ'} and @{$colN{'low_GQ'}} == 1 ) {
			$ta[ $colN{'low_GQ'}[0] ] eq 'NA' and $ta[ $colN{'low_GQ'}[0] ] = 0; 
			$ta[ $colN{'low_GQ'}[0] ] >= $opts{'minGQ'} or next SITE; 
		} elsif ( defined $colN{'low_GQ'} and @{$colN{'low_GQ'}} != 1 ) {
			; 
		} else {
			&tsmsg("[Wrn] Failed to find low_GQ for site @ta[ $colN{'CHROM'}[0], $colN{'POS'}[0] ]\n"); 
		}
	}
	print STDOUT join("\t", @ta[ $colN{'CHROM'}[0], $colN{'POS'}[0] ], @curr{qw/ref_AD_H alt_AD_H ref_AD_L alt_AD_L/})."\n"; 
}




