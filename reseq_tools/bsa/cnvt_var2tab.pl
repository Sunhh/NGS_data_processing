#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"highID:s", "lowID:s", 
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

my %colN; 
SITE: 
while (<>) {
	chomp; 
	my @ta = &splitL("\t", $_); 
	if ( $. == 1 ) {
		for (my $i=0; $i<@ta; $i++) {
			if      ($ta[$i] =~ m!^(\S+)\.AD$!) {
				if      ($1 eq $opts{'highID'}) {
					$colN{'high_AD'} = $i; 
				} elsif ($1 eq $opts{'lowID'}) {
					$colN{'low_AD'}  = $i; 
				}
			} elsif ($ta[$i] =~ m!^(\S+)\.DP$!) {
				if      ($1 eq $opts{'highID'}) {
					$colN{'high_DP'} = $i; 
				} elsif ($1 eq $opts{'lowID'}) {
					$colN{'low_DP'}  = $i; 
				}
			} elsif ($ta[$i] =~ m!^(\S+)\.GQ$!) {
				if      ($1 eq $opts{'highID'}) {
					$colN{'high_GQ'} = $i; 
				} elsif ($1 eq $opts{'lowID'}) {
					$colN{'low_GQ'}  = $i; 
				}
			} elsif ($ta[$i] =~ m!^(\S+)\.PL$!) {
				if      ($1 eq $opts{'highID'}) {
					$colN{'high_PL'} = $i; 
				} elsif ($1 eq $opts{'lowID'}) {
					$colN{'low_PL'}  = $i; 
				}
			} elsif ($ta[$i] =~ m!^CHROM$!) {
				$colN{'CHROM'} = $i; 
			} elsif ($ta[$i] =~ m!^POS$!) {
				$colN{'POS'} = $i; 
			} else {
				&tsmsg("[Wrn] Skip col-$i of [$ta[$i]]\n"); 
			}
		}
		for my $k1 (qw/high_AD high_DP low_AD low_DP CHROM POS/) {
			defined $colN{$k1} or &stopErr("[Err] failed to find column for [$k1]\n"); 
		}
		print STDOUT join("\t", qw/CHROM POS n2_RefAD_H n4_AltAD_H n1_RefAD_L n3_AltAD_L/)."\n"; 
		next SITE; 
	}
	$ta[ $colN{'high_AD'} ] =~ m!^0(,0)+$!i and next SITE; 
	$ta[ $colN{'low_AD'} ]  =~ m!^0(,0)+$!i and next SITE; 
	if ($opts{'onlyBiAllele'}) {
		$ta[ $colN{'high_AD'} ] =~ m!^\d+,\d+$! or next SITE; 
		$ta[ $colN{'low_AD'}  ] =~ m!^\d+,\d+$! or next SITE; 
	}
	if ( $ta[ $colN{'high_DP'} ] eq 'NA' ) {
		print STDERR join("\t", qw/high_AD high_DP low_AD low_DP CHROM POS high_GQ low_GQ/)."\n"; 
		print STDERR join("\t", @ta[ @colN{qw/high_AD high_DP low_AD low_DP CHROM POS high_GQ low_GQ/} ])."\n"; 
		die "$_\n"; 
	}
	my $ref_AD_H = (split(/,/, $ta[ $colN{'high_AD'} ]))[0]; $ref_AD_H =~ s!^\s+|\s+$!!g; 
	my $alt_AD_H = $ta[ $colN{'high_DP'} ]-$ref_AD_H; 
	my $ref_AD_L = (split(/,/, $ta[ $colN{'low_AD'} ]))[0];  $ref_AD_L =~ s!^\s+|\s+$!!g; 
	my $alt_AD_L = $ta[ $colN{'low_DP'} ]-$ref_AD_L; 
	my $totalDep = $ref_AD_H+$alt_AD_H+$ref_AD_L+$alt_AD_L; 
	$totalDep > 0 or next SITE; 
	$opts{'minTotalDepth'} > 0 and $totalDep < $opts{'minTotalDepth'} and next SITE; 
	$opts{'maxTotalDepth'} > 0 and $totalDep > $opts{'maxTotalDepth'} and next SITE; 
	my $refAF = ($ref_AD_H+$ref_AD_L) / $totalDep; 
	$opts{'minRefAF'} > 0 and $refAF < $opts{'minRefAF'} and next SITE; 
	$opts{'maxRefAF'} > 0 and $refAF > $opts{'maxRefAF'} and next SITE; 
	if ($opts{'minSampleDepth'} > 0) {
		$ref_AD_H+$alt_AD_H >= $opts{'minSampleDepth'} or next SITE; 
		$ref_AD_L+$alt_AD_L >= $opts{'minSampleDepth'} or next SITE; 
	}
	if ($opts{'minGQ'} > -10) {
		if (defined $colN{'high_GQ'}) {
			$ta[ $colN{'high_GQ'} ] eq 'NA' and $ta[ $colN{'high_GQ'} ] = 0; 
			$ta[$colN{'high_GQ'}] >= $opts{'minGQ'} or next SITE; 
		} else {
			&tsmsg("[Wrn] Failed to find high_GQ for site @ta[ $colN{'CHROM'}, $colN{'POS'} ]\n"); 
		}
		if (defined $colN{'low_GQ'}) {
			$ta[ $colN{'low_GQ'} ] eq 'NA' and $ta[ $colN{'low_GQ'} ] = 0; 
			$ta[$colN{'low_GQ'}] >= $opts{'minGQ'} or next SITE; 
		} else {
			&tsmsg("[Wrn] Failed to find low_GQ for site @ta[ $colN{'CHROM'}, $colN{'POS'} ]\n"); 
		}
	}
	print STDOUT join("\t", @ta[@colN{qw/CHROM POS/}], $ref_AD_H, $alt_AD_H, $ref_AD_L, $alt_AD_L)."\n"; 
}




