#!/usr/bin/perl
use strict; 
use warnings; 

!@ARGV and -t and die "perl $0 tree_files\n"; 

print <<STDOUT;
#NEXUS
Begin taxa;
STDOUT

my (%tax, @tree_list); 
my %tagH; 
while (<>) {
	chomp; 
	m/^\s*(#|$)/ and next; 
	my @ta = split(/\t/, $_); 
	my $old_fname = $ta[0]; 
	my $tag_fname = $old_fname; 
	$tag_fname =~ s!^\S+/!!; 
	$tag_fname =~ s!\.txt$!!; 
	$tag_fname =~ s!_tree$!!; 
	defined $ta[1] and $ta[1] ne '' and $tag_fname = $ta[1]; 
	defined $tagH{$tag_fname} and die "repeated tag_name [$tag_fname]\n"; 
	$tagH{$tag_fname} = 1; 
	
	open F,'<',"$old_fname" or die; 
	while (my $l = <F>) {
		chomp($l); 
		while ($l =~ s!([^\s():,]+?)_[^\s():,]+!$1!) {
			$tax{$1} //= 1; 
		}
		while ($l =~ s!\)\s*[\d\.]+!)!) {
			; 
		}
		push(@tree_list, [ $tag_fname, $l ]); 
	}
	close F; 
}
my $ntax = scalar(keys %tax); 
print <<STDOUT;
	Dimensions ntax=$ntax;
		taxlabels
STDOUT
for my $t (sort keys %tax) {
	print STDOUT "\t\t\t'$t'\n"; 
}
print STDOUT "\t\t\t;\n"; 
print STDOUT "End;\n"; 
print STDOUT "\n"; 

print STDOUT "Begin trees;\n"; 
for my $t (@tree_list) {
	print STDOUT "tree $t->[0] = $t->[1]\n"; 
}
print STDOUT "End;\n"; 

