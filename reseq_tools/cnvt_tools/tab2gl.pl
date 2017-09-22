#!/usr/bin/perl
## http://popgen.dk/angsd/index.php/Genotype_Likelihoods#Output_genotype_likelihoods
## Convert vcfTab to 'beagle genotype likelihood format' . 
use strict; 
use warnings; 
use fileSunhh; 

my $o1f = 'ind2realID_list'; 
!@ARGV and die "perl $0 in.vcf.tab > ngsadmix_input\nOut $o1f\nConvert vcfTab to 'beagle genotype likelihood format'(gl)\n"; 

my $fn = shift; 

my %geno2num = qw(
A 0
C 1
G 2
T 3
);
my $max_geno2num = 3;  

my $fh = &openFH( $fn, '<' ); 

my $header = <$fh>; 
chomp($header); 
my @hh = &splitL("\t", $header); 
print join("\t", qw/marker allele1 allele2/); 
unless ( -e $o1f ) {
	open O1,'>',"$o1f" or die; 
	for (my $i=3; $i<@hh; $i++) {
		my $iID = $i-3; 
		print "\tInd$iID\tInd$iID\tInd$iID"; 
		print O1 "Ind$iID\t$hh[$i]\n"; 
	}
	close O1; 
} else {
	for (my $i=3; $i<@hh; $i++) {
		my $iID = $i-3; 
		print "\tInd$iID\tInd$iID\tInd$iID"; 
		print STDERR "Ind$iID\t$hh[$i]\n"; 
	}
}
print "\n"; 
while (<$fh>) {
	chomp; 
	my @ta = &splitL("\t", $_); 
	my $alN = 0; 
	my $a1 = $ta[2]; 
	$a1 = uc($ta[2]); 
	($a1 ne 'N' and $a1 ne '.') and $alN ++; 
	my $a2; 
	for (my $i=3; $i<@hh; $i++) {
		$ta[$i] = uc($ta[$i]); 
		$ta[$i] =~ m!^([ATGCN.*]+)/([ATGCN.*]+)$! or die "$_\n$ta[$i]\n"; 
		my ($t1, $t2) = ($1, $2); 
		($t1 eq 'N' or $t1 eq '.') and next; 
		($t2 eq 'N' or $t2 eq '.') and next; 
		if ($a1 eq 'N' or $a1 eq '.') {
			$a1 = $t1; 
			$alN ++; 
		}
		if ( $t1 ne $a1 ) {
			unless (defined $a2) {
				$a2 = $t1; 
				$alN ++; 
			}
			$a2 ne $t1 and do { $alN ++; last; }; 
		}
		if ( $t2 ne $a1 ) {
			unless (defined $a2) {
				$a2 = $t2; 
				$alN++; 
			}
			$a2 ne $t2 and do { $alN ++; last; }; 
		}
	}
	$alN == 2 or next; 
	defined $geno2num{$a1} or do { $max_geno2num ++; $geno2num{$a1} = $max_geno2num; }; 
	defined $geno2num{$a2} or do { $max_geno2num ++; $geno2num{$a2} = $max_geno2num; }; 
	# print join("\t", "$ta[0]_$ta[1]", $a1, $a2); 
	print join("\t", "$ta[0]_$ta[1]", $geno2num{$a1}, $geno2num{$a2}); 

	for (my $i=3; $i<@hh; $i++) {
		$ta[$i] =~ m!^([ATGCN.*]+)/([ATGCN.*]+)$! or die "$_\n$ta[$i]\n"; 
		my ($t1, $t2) = ($1, $2);
		if ( $t1 =~ m/^[.N]$/ and $t2 =~ m/^[.N]$/ ) {
			print "\t0.333\t0.333\t0.333"; 
		} elsif ( $t1 eq $a1 and $t2 eq $a1 ) {
			print "\t0.999\t0.001\t0.000"; 
		} elsif ( $t1 eq $a2 and $t2 eq $a2 ) {
			print "\t0.000\t0.001\t0.999"; 
		} elsif ( ($t1 eq $a1 and $t2 eq $a2) or ($t1 eq $a2 and $t2 eq $a1) ) {
			print "\t0.000\t0.999\t0.000"; 
		} else {
			die "i=$i $ta[$i] $t1 $t2 [$a1 $a2]\n"; 
		}
	}
	print "\n"; 
}
close($fh); 

