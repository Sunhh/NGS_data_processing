#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use SNP_tbl; 

!@ARGV and die "perl $0 in.cols > out.vcf\n"; 

my $st = SNP_tbl->new(); 

my $is_noSort = 0; 

my %dblist; 
my %bb2b; 
{
my @aa = (
[qw/W A T/], 
[qw/S C G/],
[qw/M A C/],
[qw/K G T/],
[qw/R A G/],
[qw/Y C T/], 
[qw/A A A/], 
[qw/T T T/], 
[qw/G G G/], 
[qw/C C C/]
); 
for my $tr (@aa) {
	my @bb = @$tr; 
	$dblist{$bb[0]} = [$bb[1], $bb[2]]; 
	$bb2b{"$bb[1]$bb[2]"} = $bb[0]; 
	$bb2b{"$bb[2]$bb[1]"} = $bb[0]; 
}
}


# .ped file format : 
#   pedID indvID FatherID MotherID sexTag caseTag genotypes(A=1, C=2, G=3, T=4, Miss=0)
my %allele2num = qw(
A 1
C 2
G 3
T 4
N 0
); 

print STDOUT <<HH;
##fileformat=VCFv4.1
##source=SNP_tbl
##FILTER=<ID=PASS,Description="Passed variant FILTERs">
HH

my $snpCol = 3; 
my (@header, @data, @loci); 
while (<>) {
	$. % 100e3 == 1 and &tsmsg("[Msg] $. lines.\n"); 
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($. == 1) {
		@header = @ta; 
		print STDOUT join("\t", '#CHROM', qw/POS ID REF ALT QUAL FILTER INFO FORMAT/, @header[$snpCol .. $#header])."\n"; 
		next; 
	}

	my $chr_ori = $ta[0]; 
	my $chr = $chr_ori; 
	$chr =~ s!^(WM97(?i:v\d*)?_)?Chr0*!!i; 
	$chr eq '' and $chr = 20; 
	$chr =~ m/^\d+$/ or $chr = uc($chr); 
	# $chr =~ m/^\d+$/ or &stopErr("[Err] chr=$chr not a number\n"); 

	my $posi = $ta[1]; 

	my @allele_arr; 
	my ($refBase, %baseNum, $cntNum); 

	$refBase = $ta[2]; 
	if ( length($refBase) > 1 ) {
		if ($refBase =~ m/^([ATGC]+)$/) {
			my @ta = sort { $a cmp $b } split(//, $refBase); 
			$refBase = $ta[0]; 
#		} else {
#			defined $bb2b{$refBase} or die "refBase [$refBase]\n"; 
#			$refBase = $bb2b{$refBase}; 
#			defined $dblist{$refBase} or die "db refBase [$refBase]\n"; 
#			$refBase = $dblist{$refBase}[0]; 
		}
	}

	$cntNum = 1; 
	for my $bb ( @ta[$snpCol .. $#ta] ) {
		$bb = uc($bb); 
		if ( $bb =~ m!^([ATGC]{1,2})$! ) {
			$bb = $st->SingleChar($bb); 
			defined $dblist{ $bb } or do { push( @allele_arr, ['.', '.'] ); next; }; 
			push( @allele_arr, [ @{$dblist{ $bb }} ] ); 
		} elsif ( $bb =~ m!^([ATGC])\+([ATGCN]+)$! ) {
			$bb = "$1$2"; 
			push( @allele_arr, [ $bb, $bb ] ); 
		} elsif ( $bb eq '*' ) {
			push( @allele_arr, [ '*', '*' ] ); 
		} elsif ( $bb =~ m!^([ATGC])(\*)$! ) {
			push( @allele_arr, [ $1, $2 ] ); 
		} elsif ( $bb =~ m!^(\*)([ATGC])$! ) {
			push( @allele_arr, [ $2, $1 ]); 
		} elsif ( $bb =~ m!^([ATGC])([ATGC])\+([ATGCN]+)$! ) {
			push( @allele_arr, [ $1, "$2$3" ] ); 
		} elsif ( $bb =~ m!^\+([ATGCN]+)$! ) {
			push( @allele_arr, [ $1, $1 ] ); 
		} else {
			push( @allele_arr, [ '.', '.' ] ); 
		}
		for my $b2 ( @{$allele_arr[-1]} ) {
			$b2 eq '.' and next; 
			defined $baseNum{$b2} and next; 
			# defined $refBase or $refBase = $b2; 
			$baseNum{$b2} = $cntNum; 
			$cntNum ++; 
		}
	}
	defined $refBase or die "$_\n"; 
	$cntNum >= 2 or die "2: $_\n"; 
	$baseNum{$refBase} = 0; 
	my @allele_base = sort { $baseNum{$a} <=> $baseNum{$b} } keys %baseNum; 
	@allele_base == 1 and $allele_base[1] = '.'; 
	for (my $i=0; $i<@allele_base; $i++) {
		$allele_base[$i] eq '.' and next; 
		$baseNum{$allele_base[$i]} = $i; 
	}

	my @allele_str; 
	for my $tr1 (@allele_arr) {
		if ( $tr1->[0] eq '.' ) {
			push( @allele_str, './.' ); 
		} else {
			$is_noSort or ( $baseNum{$tr1->[0]} > $baseNum{$tr1->[1]} and ( $baseNum{$tr1->[0]}, $baseNum{$tr1->[1]} ) = ( $baseNum{$tr1->[1]}, $baseNum{$tr1->[0]} ) ); 
			push( @allele_str, join('/', $baseNum{$tr1->[0]}, $baseNum{$tr1->[1]}) );
		}
	}

	# print STDOUT join( "\t", $chr_ori, $posi, "s${chr}_${posi}", $allele_base[0], join(',', @allele_base[1 .. $#allele_base]), qw/. PASS . GT/, @allele_str)."\n"; 
	print STDOUT join( "\t", $chr_ori, $posi, "s${chr}_${posi}", $refBase, join(',', @allele_base[1 .. $#allele_base]), qw/. PASS . GT/, @allele_str)."\n"; 

}

sub geno2num {
	my $tb = shift; 
	$tb = $st->SingleChar($tb); 
	if (defined $allele2num{$tb}) {
		return ($allele2num{$tb}, $allele2num{$tb}); 
	} elsif ( defined $dblist{$tb} ) {
		my @b1 = @{$dblist{ $tb }}; 
		return ($allele2num{$b1[0]}, $allele2num{$b1[1]}); 
	} else {
		return (0, 0); 
	}
}

