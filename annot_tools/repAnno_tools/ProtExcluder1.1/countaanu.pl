#! /usr/bin/perl -w
use strict; 
use warnings; 

my $usage = "countaanu.pl blastx.mt file\n";

# to calculate the number of each amino acids among matched ones

my $AAstring='ABCDEFGHIKLMNPQRSTVWYZ'; 
my %a2num; 
{
	my @ta = split(//, $AAstring); 
	for (my $i=0; $i<@ta; $i++) {
		$a2num{$ta[$i]} = $i; 
	}
}

if (@ARGV < 1) {die $usage;}

open(FA, "$ARGV[0]") || die $usage;

my $seq = "";
my $ct; 
my %nu; 
while (<FA>) {
	if (/>\s*(.+)/) {
		my $hh = $1; 
		if ($seq ne '') {
			my @sym = split(//, $seq);
			$ct = 0;
			foreach my $sym (@sym) {
				if ( defined $a2num{$sym} ) {
					$nu{ $a2num{$sym} } ++; 
					$ct ++; 
				}
			}
			if ($ct > 0) {
				printf "%05d ",$ct;
			}
			&report(\%nu);
		}
		my $i = 0;
		while ($i < 22) {
			$nu{$i} = 0;
			$i ++;
		}
		printf "\n>%s ", $hh;
		$seq = "";
	} else {
		chomp;
		$seq .= $_;
	}
}
close FA;


if ($seq ne '') {
	my @sym = split(//, $seq);
	$ct = 0;
	foreach my $sym (@sym) {
		if ( defined $a2num{$sym} ) {
			$nu{ $a2num{$sym} } ++; 
			$ct ++; 
		}
	}
	if ($ct > 0) {
		printf "%05d ",$ct;
	}
	&report(\%nu);
}
print "\n";

sub report {
	my ($rh) = @_; 
	for my $key (sort { $rh->{$b} <=> $rh->{$a} || $a cmp $b } keys %$rh) {
		printf("%02d ", $rh->{$key});
	}
}
