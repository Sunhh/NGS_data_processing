#!/usr/bin/perl
use strict; 
use warnings; 
@ARGV >= 2 or die "perl $0 raw.xmap new.xmap\n"; 
my $fn = shift; 
my $new_fn = shift; 
$fn =~ s/\.xmap$//; 
$new_fn =~ s/\.xmap$//; 

open F,'<',"${fn}.xmap" or die; 
open O,'>',"${new_fn}.xmap" or die; 
while (<F>) {
	if (m/^#/) {
		s/^(# Reference Maps From:\s*)\S+\.cmap/${1}${new_fn}_r.cmap/; 
		s/^(# Query Maps From:\s*)\S+\.cmap/${1}${new_fn}_q.cmap/; 
		print O; 
	} else {
		print O; 
	}
}
close F; 
system "rm ${fn}.xmap"; 
system "mv ${fn}_q.cmap ${new_fn}_q.cmap"; 
system "mv ${fn}_r.cmap ${new_fn}_r.cmap"; 

