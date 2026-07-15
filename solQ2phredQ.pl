#!/usr/bin/perl
use warnings;
use strict;

-t and !@ARGV and die "perl $0 in_sol.fq\n"; 

# rawS=|>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|
# newS=|# !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]

my $tn=0;
while (my $rdID=<>) {
	$tn % 1000000 == 1 and &tmsg("$tn lines.\n"); 
	$rdID =~ m/^@/ or die "Format wrong:$_\n";
	my $rdSeq=<>; <>;
	my $rdQua=<>;
	chomp($rdQua);
	$rdQua =~ tr!>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[ \\ ]^_`abcdefghijklmnopqrstuvwxyz \{ |}~!#####$%&'()*+,-./0123456789:;< \= >?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[ \\ ]^_!; 

	print STDOUT $rdID . ${rdSeq} . "+\n" . $rdQua . "\n";
	$tn ++;
}
sub tmsg {
	my $tt=scalar(localtime()); 
	warn join("", "$tt ", @_); 
}

