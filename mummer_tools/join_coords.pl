#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"maxDist1:i", 
	"maxDist2:i", 
	"help!", 
); 
my $help_txt = <<HH; 

perl $0 pb_ngs.coords > pb_ngs.coords.jn

-maxDist1    [40e3] Distance to join two neighboring blocks. 
-maxDist2    [40e3] Distance to join two neighboring blocks.
HH

$opts{'maxDist1'} //= 40e3; 
$opts{'maxDist2'} //= 40e3; 

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

# /Data/Sunhh/compare_NGS_pb/05_alignByMum/WM97pbV0.ctg.fa /Data/Sunhh/compare_NGS_pb/05_alignByMum/WM97_v1.scf.fa
# NUCMER
# 
#     [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]
# ===============================================================================================================================
#    16657   113511  |        1    95242  |    96855    95242  |    73.81  |  8241642    96975  |     1.18    98.21  | WM97pbV0_000000F   WM97_scaffold1452
#   128093   129621  |        1     1552  |     1529     1552  |    96.14  |  8241642    10579  |     0.02    14.67  | WM97pbV0_000000F   WM97_scaffold16934
#   134474   134996  |      524        1  |      523      524  |    99.81  |  8241642      524  |     0.01   100.00  | WM97pbV0_000000F   WM97_scaffold10350



my (%blk_F, %blk_R, %len1, %len2); 
while (<>) {
	$. == 1 and $_ =~ m!^(/|\.)! and next; 
	$_ =~ m!^NUCMER! and next; 
	m!^\s*$! and next; 
	chomp; 
	s!^\s+!!g; 
	my @ta = split(/\s+/, $_); 
	$ta[0] eq '[S1]' and next; 
	$ta[0] =~ m!^\=! and next; 
	if ( $ta[3] < $ta[4] ) {
		push(@{$blk_F{$ta[17]}{$ta[18]}}, [ 
				@ta[0,1,3,4],                                                           # S1,E1, S2,E2
				int(($ta[1]-$ta[0]+1)*$ta[9]/100), int(($ta[4]-$ta[3]+1)*$ta[9]/100),   # Match1,Match2
				[ [@ta[0,1,3,4]] ],                                                     # SE_sets; [S1,E1,S2,E2]...
				$ta[1]-$ta[0]+1, $ta[4]-$ta[3]+1                                        # BlkLen1,BlkLen2
			]
		); 
	} elsif ( $ta[3] > $ta[4] ) {
		push(@{$blk_R{$ta[17]}{$ta[18]}}, [ 
				(@ta[0,1], -$ta[3], -$ta[4]), 
				int(($ta[1]-$ta[0]+1)*$ta[9]/100), int(($ta[3]-$ta[4]+1)*$ta[9]/100),
				[ [ @ta[0,1], -$ta[3], -$ta[4] ] ], 
				$ta[1]-$ta[0]+1, $ta[3]-$ta[4]+1
		       ]); 
	} else {
		die "bad line: $_\n"; 
	}
	$len1{$ta[17]} //= $ta[11]; 
	$len2{$ta[18]} //= $ta[12]; 
}

for my $k1 (sort {$len1{$b} <=> $len1{$a} || $a cmp $b } keys %len1) {
	for my $k2 (sort { $len2{$b} <=> $len2{$a} || $a cmp $b } keys %{$blk_F{$k1}}) {
		&merge_blks($blk_F{$k1}{$k2}, $opts{'maxDist1'}, $opts{'maxDist2'}); 
		for my $t1 (@{$blk_F{$k1}{$k2}}) {
			my @ta1 = @$t1; 
			print STDOUT join("\t", @ta1[0..3], 
				sprintf("%.2f", $ta1[4]/($ta1[1]-$ta1[0]+1)*100), 
				sprintf("%.2f", $ta1[5]/($ta1[3]-$ta1[2]+1)*100), 
				sprintf("%.2f", $ta1[4]/$ta1[7]*100), 
				sprintf("%.2f", $ta1[5]/$ta1[8]*100), 
				$ta1[1]-$ta1[0]+1, $ta1[3]-$ta1[2]+1, 
				"+", $len1{$k1}, $len2{$k2}, $k1, $k2, 
				join(";", map { "$_->[0]-$_->[1]:$_->[2]-$_->[3]" } @{$ta1[6]})
			)."\n"; 
		}
	}
	for my $k2 (sort { $len2{$b} <=> $len2{$a} || $a cmp $b } keys %{$blk_R{$k1}}) {
		&merge_blks($blk_R{$k1}{$k2}, $opts{'maxDist1'}, $opts{'maxDist2'}); 
		for my $t1 (@{$blk_R{$k1}{$k2}}) {
			my @ta1 = @$t1; 
			print STDOUT join("\t", ( @ta1[0,1], -$ta1[3], -$ta1[2] ), 
				sprintf("%.2f", $ta1[4]/($ta1[1]-$ta1[0]+1)*100), 
				sprintf("%.2f", $ta1[5]/($ta1[3]-$ta1[2]+1)*100), 
				sprintf("%.2f", $ta1[4]/$ta1[7]*100), 
				sprintf("%.2f", $ta1[5]/$ta1[8]*100), 
				$ta1[1]-$ta1[0]+1, $ta1[3]-$ta1[2]+1, 
				"-", $len1{$k1}, $len2{$k2}, $k1, $k2, 
				join(";", map { my @bb=@$_; $bb[2]*=-1; $bb[3]*=-1; "$bb[0]-$bb[1]:$bb[3]-$bb[2]"; } @{$ta1[6]})
			)."\n"; 
		}
	}
}

sub merge_blks {
	my ($blkR, $maxDist1, $maxDist2) = @_; 
	@$blkR = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3] } @$blkR; 
	while ( &shortenBlk($blkR, $maxDist1, $maxDist2) ) {
		1; 
	}
	return($blkR); 
}# sub merge_blks() 

sub shortenBlk {
	my ($blkR, $maxDist1, $maxDist2) = @_; 
	my $changed = 0; 
	for (my $i=1; $i<@$blkR; $i++) {
		my $j= $i-1; 
		$blkR->[$i][0] > $blkR->[$j][1] or next; 
		$blkR->[$i][2] > $blkR->[$j][3] or next; 
		$blkR->[$i][0]-$blkR->[$j][1]-1 <= $maxDist1 or next; 
		$blkR->[$i][2]-$blkR->[$j][3]-1 <= $maxDist2 or next; 
		$blkR->[$j] = [ $blkR->[$j][0], $blkR->[$i][1], $blkR->[$j][2], $blkR->[$i][3], # S1,E1, S2,E2
			$blkR->[$j][4]+$blkR->[$i][4], $blkR->[$j][5]+$blkR->[$i][5],           # Match1,Match2
			[ @{$blkR->[$j][6]}, @{$blkR->[$i][6]} ],                               # SE_sets; [S1,E1,S2,E2]...
			$blkR->[$j][7]+$blkR->[$i][7], $blkR->[$j][8]+$blkR->[$i][8]
		]; 
		splice(@$blkR, $i, 1); 
		$changed = 1; 
		last; 
	}
	return($changed); 
}# sub shortenBlk() 

