#!/usr/bin/perl
use strict; 
use warnings; 
my %jn; 
my %go; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if (defined $ta[11] and $ta[11] =~ m!^IPR!) {
		push(@{$jn{$ta[0]}}, "$ta[11]:$ta[12]"); 
		push(@{$go{$ta[0]}}, ""); 
		# defined $ta[13] and $ta[13] =~ m!^GO! or die "GO!\n$_\n"; 
	}
	if ( defined $ta[13] and $ta[13] =~ m!^GO! ) {
		( defined $ta[11] and $ta[11] =~ m!^IPR! ) or die "GO?\n$_\n"; 
		$go{$ta[0]}->[-1] = $ta[13]; 
	}
}


for (keys %jn) {
	# defined $go{$_} or $go{$_} = []; 
	print STDOUT join("\t", $_, join("; ", @{$jn{$_}}), join("; ", @{$go{$_}}))."\n"; 
}

