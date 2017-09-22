#!/usr/bin/perl -w 
use strict; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"saveLines:i", 
); 

!@ARGV and die <<INFO; 

perl $0 <in_list>  <in_sites>
################################################################################
# In_list file format : 
#   RegionName        RegionChrID       RegionChrStart        RegionChrEnd
# In_sites file format : 
#   ChrID             ChrPosition       
################################################################################

-saveLines       [0]

INFO

# In_list file format : 
#   RegionName        RegionChrID       RegionChrStart        RegionChrEnd
# In_sites file format : 
#   ChrID             ChrPosition       

my $inListF = shift; 
my $inSiteF = shift; 

my $inLFH = &openFH($inListF, '<'); 
my $inSFH = &openFH($inSiteF, '<'); 

my %need_region; 
while (<$inLFH>) {
	chomp; 
	/^\s*$/ and next; 
	/^[^\S\t]*#/ and next; 
	my @ta = split(/\t/, $_); 
	$ta[2] > $ta[3] and @ta[2,3] = @ta[3,2]; 
	push(@{$need_region{$ta[1]}}, [$ta[0], $ta[2], $ta[3]]); 
}
close ($inLFH); 

# Sort needed regions. 
for my $r1 (values %need_region) {
	@$r1 = sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] || $a->[0] <=> $b->[0] || $a->[0] cmp $b->[0]} @$r1; 
}

if (defined $opts{saveLines}) {
	while (<$inSFH>) {
		print "RegionID\t" . $_; 
		$. < $opts{saveLines} or last; 
	}
}else {
	print STDOUT "RegionID\tLine\n"; 
}
while (<$inSFH>) {
	$. % 1000000 == 1 and print STDERR "[Stat] $. lines over.\n"; 
	chomp; /^\s*$/ and next; 
	my ($cid, $cloc) = split(/\t/, $_); 
	if (defined $need_region{$cid}) {
		CHK_REGION: 
		for my $r1 (@{$need_region{$cid}}) {
			my ($name, $s, $e) = @$r1; 
			if ($s > $cloc) {
				last CHK_REGION; 
			} elsif ($s <= $cloc and $e >= $cloc) {
				print STDOUT "$name\t$_\n"; 
			}
		}
	}
}
close ($inSFH); 


sub openFH {                                                                                                                                                                                                                                 
	my $f = shift;
	my $type = shift; (defined $type and $type =~ /^[<>|]+$/) or $type = '<';
	local *FH;
	open FH,$type,"$f" or die "Failed to open file [$f]:$!\n";
	return *FH;
}# end sub openFH


