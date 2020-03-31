#!/usr/bin/perl
use strict; 
use warnings; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 

# S1	E1	S2	E2	ID1_1	ID1_2	ID2_1	ID2_2	Span1	Span2	Strand	Size1	Size2	Name1	Name2	Block_Number	Block_Loci
# 16657	113511	1	95242	73.81	73.81	73.81	73.81	96855	95242	+	8241642	96975	WM97pbV0_000000F	WM97_scaffold1452	1	16657-113511:1-95242
# 128093	139235	1	10579	20.3	21.68	96.83	96.83	11143	10579	+	8241642	10579	WM97pbV0_000000F	WM97_scaffold16934	2	128093-129621:1-1552;138429-139235:9764-10579
# 134474	134996	1	524	99.81	99.81	99.81	99.81	523	524	-	8241642	524	WM97pbV0_000000F	WM97_scaffold10350	1	134474-134996:1-524


my (%blk1, %blk2); 
my (%len1, %len2); 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[0] eq 'S1' and next; 
	$ta[0] eq 'R_S' and next; 
	push(@{$blk1{$ta[13]}}, [@ta[0,1]]); 
	push(@{$blk2{$ta[14]}}, [@ta[2,3]]); 
	$len1{$ta[13]} //= $ta[11]; 
	$len2{$ta[14]} //= $ta[12]; 
}

for my $k1 (sort {$len1{$b} <=> $len1{$a}} keys %len1) {
	my $m1 = $ms_obj->mergeLocBlk( $blk1{$k1}, 'dist2join' => 1 ); 
	for my $tse (@$m1) {
		print join("\t", "ID1", $k1, @{$tse}[0,1], $tse->[1]-$tse->[0]+1)."\n"; 
	}
}

for my $k2 (sort {$len2{$b} <=> $len2{$a}} keys %len2) {
	my $m1 = $ms_obj->mergeLocBlk( $blk2{$k2}, 'dist2join' => 1 ); 
	for my $tse (@$m1) {
		print join("\t", "ID2", $k2, @{$tse}[0,1], $tse->[1]-$tse->[0]+1)."\n"; 
	}
}



