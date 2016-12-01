#!/usr/bin/perl -w
# Cannot handle insertions. 
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"noHeader!", 
); 

my $help_txt = <<HH; 

perl $0 in.vcf.tab > in.vcf.cols

# Mask indel genotypes to N, except single '*'. 
# Only accept bi-allele genotype for cols-format (single-character) genotype. 

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

my %dblist; 
{
	my @aa = (
	[qw/W A T/], 
	[qw/S C G/],
	[qw/M A C/],
	[qw/K G T/],
	[qw/R A G/],
	[qw/Y C T/]
	); 
	for my $tr (@aa) {
		my @bb = @$tr; 
		$dblist{$bb[0]} = [$bb[1], $bb[2]]; 
	}
}


unless ( defined $opts{'noHeader'} ) {
	my $l = <>; 
	print STDOUT $l; 
}

my %used; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	for my $tb ( @ta[2 .. $#ta] ) {
		$tb = uc($tb); 
		if ( $tb =~ s!^([ATGC\*])/\1$!$1! ) {
		} elsif ( $tb =~ s!^([ATGC\*])/([ATGC\*])$!$1$2! ) {
		} elsif ( $tb eq './.' or $tb eq 'N/N' ) {
			$tb = 'N'; 
		} elsif ( $tb =~ m!^[ATGC\*]+/[ATGC\*]+$! ) {
			unless ( defined $used{'skip_geno'}{$tb} ) {
				&tsmsg("[Wrn] Mask genotype [$tb]\n"); 
				$used{'skip_geno'}{$tb} = 1; 
			}
			$tb = 'N'; 
		} elsif ( $tb =~ m!^[ATGCN\*]$! ) {
			; 
		} elsif ( defined $dblist{$tb} ) {
			$tb = $dblist{$tb}[0] . $dblist{$tb}[1]; 
		} else {
			unless ( defined $used{'skip_geno'}{$tb} ) {
				&tsmsg("[Wrn] Mask genotype [$tb]\n"); 
				$used{'skip_geno'}{$tb} = 1; 
			}
			$tb = 'N'; 
		}
	}
	print STDOUT join("\t", @ta)."\n"; 
}

