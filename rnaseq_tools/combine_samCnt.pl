#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

-t and !@ARGV and die "perl $0 file_list_to_combine\n"; 

my $list = shift; 

my @out_txt; 
my %fh; 
$fh{'F'} = &openFH( $list, '<' ); 
&tsmsg("[Msg] Reading list [$list]\n"); 
my @fn_lis; 
while (&wantLineC($fh{'F'})) {
	my @ta = &splitL("\t", $_); 
	push(@fn_lis, $ta[0]); 
}
close($fh{'F'}); 

$fh{'C1'} = &openFH( $fn_lis[0], '<' ); 
&tsmsg("[Msg] Reading file [$fn_lis[0]]\n"); 
while (&wantLineC($fh{'C1'})) {
	my @ta = &splitL("\t", $_); 
	push(@out_txt, [ $ta[0], @ta[1 .. $#ta] ]); 
}
close $fh{'C1'}; 

for ( my $i=1; $i<@fn_lis; $i++ ) {
	my $fh = &openFH( $fn_lis[$i], '<' ); 
	&tsmsg("[Msg] Reading file [$fn_lis[$i]]\n"); 
	my $ln = -1; 
	while ( &wantLineC($fh) ) {
		$ln ++; # In this way, line number only corresponds to good lines. 
		my @ta = &splitL("\t", $_); 
		$#{$out_txt[$ln]} == $i or &stopErr("[Err] Not equal cols at @ta\n"); 
		$out_txt[$ln][0] eq $ta[0] or &stopErr("[Err] Different first column at @ta\n"); 
		push( @{$out_txt[$ln]}, @ta[1 .. $#ta] ); 
	}
	close($fh); 
}

&tsmsg("[Msg] Output lines.\n"); 
for (@out_txt) {
	print STDOUT join("\t", @$_)."\n"; 
}

&tsmsg("[Rec] $0 done.\n"); 

