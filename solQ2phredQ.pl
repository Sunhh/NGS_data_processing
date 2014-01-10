#!/usr/bin/perl
use warnings;
use strict;

-t and !@ARGV and die "perl $0 in_sol.fq\n"; 

# rawS=|>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|
# newS=|# !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]

#my %h;
#my @rawS = (); 
#my @newS = (); 
#for my $qual_num (-2 .. 70) {
#	my $acsii_sol = chr( $qual_num + 64 ); 
#	my $acsii_san = chr( $qual_num + 33 ); 
#	$qual_num <= 0 and $acsii_san = chr( 2 + 33 ); 
#	$h{$acsii_sol} = $acsii_san;
#	push(@rawS, $acsii_sol); 
#	push(@newS, $acsii_san); 
#}
#print "rawS=|",join("",@rawS),"\nnewS=|",join("",@newS),"\n"; 
#exit 1;

my $tn=0;
while (my $rdID=<>) {
	$tn % 1000000 == 1 and &tmsg("$tn lines.\n"); 
	$rdID =~ m/^@/ or die "Format wrong:$_\n";
	my $rdSeq=<>; <>;
	my $rdQua=<>;
	chomp($rdQua);
#	my $newQua = &sol2san(\$rdQua);
	$rdQua =~ tr!>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[ \\ ]^_`abcdefghijklmnopqrstuvwxyz \{ |}~!#####$%&'()*+,-./0123456789:;< \= >?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[ \\ ]^_!; 

	print STDOUT $rdID . ${rdSeq} . "+\n" . $rdQua . "\n";
	$tn ++;
}
sub tmsg {
	my $tt=scalar(localtime()); 
	warn join("", "$tt ", @_); 
}

#sub sol2san {
#	my $qR=shift(@_);
#	my $back = '';
#	while ($$qR =~ /(.)/g) {
#		$back .= $h{$1};
#	}
#	return $back;
#}

