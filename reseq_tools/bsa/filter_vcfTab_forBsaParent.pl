#!/usr/bin/perl
use strict; 
use warnings; 
use Getopt::Long; 
use LogInforSunhh; 
my %opts; 
GetOptions(\%opts, 
	"col_noN:i@", # 
	"col_homo:i@", # 
	"col_diffPair:s@", 
	"col_samePair:s@", 
	"col_diff2Ref:s@", 
	"help!", 
	"noHeader!", 
); 

-t and !@ARGV and die "perl $0 input.vcfTab -col_noN 3 -col_noN 4 -col_homo 3 -col_homo 4 -col_diffPair 3,4 -col_diffPair 3,5 -col_samePair 4,5 > input.vcfTab.filtered\n"; 

unless ($opts{'noHeader'}) {
	my $hh = <>; 
	print $hh; 
}

my ($chk_noN, $chk_homo, $chk_diffP, $chk_sameP, $chk_diff2Ref) = (0,0,0,0,0); 
my (@cn_noN, @cn_homo, @cn_diffP, @cn_sameP, @cn_diff2Ref); 
if (defined $opts{'col_noN'}) {
	push(@cn_noN, @{$opts{'col_noN'}}); 
	$chk_noN = 1; 
}
if (defined $opts{'col_homo'}) {
	push(@cn_homo, @{$opts{'col_homo'}}); 
	$chk_homo = 1; 
}
if (defined $opts{'col_diff2Ref'}) {
	push(@cn_diff2Ref, @{$opts{'col_diff2Ref'}}); 
	$chk_diff2Ref = 1; 
}
if (defined $opts{'col_diffPair'}) {
	for my $t1 (@{$opts{'col_diffPair'}}) {
		$t1 =~ m!^(\d+),(\d+)$! or die "bad input [$t1]\n"; 
		push(@cn_diffP, [$1, $2]); 
	}
	$chk_diffP = 1; 
}
if (defined $opts{'col_samePair'}) {
	for my $t1 (@{$opts{'col_samePair'}}) {
		$t1 =~ m!^(\d+),(\d+)$! or die "bad input [$t1]\n"; 
		push(@cn_sameP, [$1, $2]); 
	}
	$chk_sameP = 1; 
}

LINE: 
while (<>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	if ($chk_noN == 1) {
		for my $t2 (@ta[@cn_noN]) {
			$t2 =~ m!(^\./\.$|N|n)! and next LINE; 
		}
	}
	if ($chk_homo == 1) {
		for my $t2 (@ta[@cn_homo]) {
			$t2 =~ m!^[ATGCatgc]+$! or $t2 =~ m!^([ATGCatgc]+)/\1$! or next LINE; 
		}
	}
	if ($chk_diff2Ref == 1) {
		for my $t2 (@ta[@cn_diff2Ref]) {
			$t2 eq "$ta[2]/$ta[2]" and next LINE; 
		}
	}
	if ($chk_diffP == 1) {
		for my $t3 (@cn_diffP) {
			my ($cn1, $cn2) = @$t3; 
			$ta[$cn1] eq $ta[$cn2] and next LINE; 
		}
	}
	if ($chk_sameP == 1) {
		for my $t3 (@cn_sameP) {
			my ($cn1, $cn2) = @$t3; 
			$ta[$cn1] eq $ta[$cn2] or next LINE; 
		}
	}

}


