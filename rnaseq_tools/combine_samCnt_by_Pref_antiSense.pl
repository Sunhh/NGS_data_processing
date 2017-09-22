#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

!@ARGV and die "perl $0 prefix_of_sample [./] [P1|P3]\n"; 

my $pref = shift; 
my $dir = shift; 
my $pTag = shift; 
$dir //= './'; 
$pTag //= 'P1'; 
$pTag =~ m/^(P1|P3)$/ or die "pTag=[$pTag]\n"; 


my @suff_arr; 
if ($pTag eq 'P1') {
	push( @suff_arr, ["P1g_allMisR0p06", "_toP1.allMisR0p06.anti.cnt"] ) ; 
	push( @suff_arr, ["P1g_P1spec",      "_toP1.specMisR0p01_inP1.anti.cnt"] ) ; 
	push( @suff_arr, ["P1g_P3spec",      "_toP3.specMisR0p01_inP1.anti.cnt"] ) ; 
	push( @suff_arr, ["P1g_cmmn",        "_toP1.cmmnMisR0p06_inP1.anti.cnt"] ) ; 
	# Combine subsets. 
	#  "P1g_total", 'P1g_P1spec' + 'P1g_P3spec' + 'P1g_cmmn'
} elsif ($pTag eq 'P3') {
	push( @suff_arr, ["P3g_allMisR0p06", "_toP3.allMisR0p06.anti.cnt"] ) ; 
	push( @suff_arr, ["P3g_P1spec",      "_toP1.specMisR0p01_inP3.anti.cnt"] ) ; 
	push( @suff_arr, ["P3g_P3spec",      "_toP3.specMisR0p01_inP3.anti.cnt"] ) ; 
	push( @suff_arr, ["P3g_cmmn",        "_toP3.cmmnMisR0p06_inP3.anti.cnt"] ) ; 
	#  "P3g_total", 'P3g_P1spec' + 'P3g_P3spec' + 'P3g_cmmn'
} else {
	die "fasdf\n"; 
}

my @out_txt; 
my @out_header; 
@out_header = ("Gene_ID"); 
my @fn_lis; 
for my $suff (@suff_arr) {
	push(@fn_lis, $dir . "${pref}" . $suff->[1]); 
	push(@out_header, "$suff->[0]_$pref"); 
}

my %fh; 

$fh{'C1'} = &openFH( $fn_lis[0], '<' ); 
&tsmsg("[Msg] Reading file [$fn_lis[0]]\n"); 
while (&wantLineC($fh{'C1'})) {
	my @ta = &splitL("\t", $_); 
	# push(@out_txt, [ $ta[0], @ta[1 .. $#ta] ]); 
	push(@out_txt, [ $ta[0] ]); 
}
close $fh{'C1'}; 

for ( my $i=0; $i<@fn_lis; $i++ ) {
	my $fh = &openFH( $fn_lis[$i], '<' ); 
	&tsmsg("[Msg] Reading file [$fn_lis[$i]]\n"); 
	my $ln = -1; 
	while ( &wantLineC($fh) ) {
		$ln ++; # In this way, line number only corresponds to good lines. 
		my @ta = &splitL("\t", $_); 
		$#{$out_txt[$ln]} == $i or &stopErr("[Err] Not equal cols at @ta\n"); 
		$out_txt[$ln][0] eq $ta[0] or &stopErr("[Err] Different first column at @ta\n"); 
		push( @{$out_txt[$ln]}, $ta[1] ); 
	}
	close($fh); 
}

if ($pTag eq 'P1') {
	push(@out_header, "P1g_total_$pref"); 
} elsif ($pTag eq 'P3') {
	push(@out_header, "P3g_total_$pref"); 
} else {
	die "aaa\n"; 
}

for (my $i=1; $i<@out_txt; $i++) {
	$out_txt[$i][ $#{$out_txt[$i]}+1 ] = $out_txt[$i][2] + $out_txt[$i][3] + $out_txt[$i][4]; 
}

&tsmsg("[Msg] Output lines.\n"); 
print STDOUT join("\t", @out_header)."\n"; 
for (my $i=1; $i<@out_txt; $i++) {
	print STDOUT join("\t", @{$out_txt[$i]})."\n"; 
}

&tsmsg("[Rec] $0 done.\n"); 

