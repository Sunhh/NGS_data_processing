#!/usr/bin/perl
use strict; 
use warnings; 

my $pd0 = shift; 
my @all_sub = &list_sub($pd0); 
for (@all_sub) {
	print "$_\n"; 
}

sub list_sub {
	# Skip links; 
	my ($pd) = @_; 
	# warn "pd=$pd\n"; 
	my @subB = ($pd); 
	if (! -d $pd ) {
		# $pd is not a folder; return $pd itself; 
		return(@subB); 
	}
	# $pd is a folder, get the children of $pd; 
	opendir D0,$pd or die "failed to opendir [$pd]\n"; 
	my @sub1 = map { "$pd/$_" } grep { $_ !~ m!^\.\.?$! } readdir(D0); 
	push(@subB, @sub1); 
	# There might be subdirs in @sub1, so check it. 
	for my $sd1 (@sub1) {
		my @sub2 = &list_sub($sd1); 
		push(@subB, @sub2[1..$#sub2]); 
	}
	return(@subB); 
}# list_sub() 


