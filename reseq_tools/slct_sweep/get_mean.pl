#!/usr/bin/perl 
# ChrID   WindS   WindE   WindL   BpCnt   pi_gAll theta_gAll      tajima_D_gAll   pi_g4   theta_g4        tajima_D_g4     pi_g3   theta_g3        tajima_D_g3     pi_g2   theta_g2        tajima_D_g2     pi_g2a3 theta_g2a3      tajima_D_g2a3
# chr1    1       100000  100000  27657   59.8128321337606        136.979297431795        -1.9207331754534        99.3359589967869        149.831476834393        -1.37327869278475       58.211999739167 76.7064104362463        -1.20362280767511    27.3923716438263        54.2375154752775        -1.92707130256499       45.6253736331524        71.6667508015005        -1.36068427374064
# chr1    10001   110000  100000  27822   76.8549472740237        151.162969364449        -1.67690170040421       126.716125250104        166.35899959551 -0.971428294573386      67.5341894517592        94.3806985551971        -1.421720604758      32.92761442038  62.8951470534908        -1.85817356513585       56.4772732872807        85.3851387437019        -1.26973636792564

use strict; 
use warnings; 

!@ARGV and die "perl $0 apple_noN_w100ks10k.grpAll_grp4_grp3_grp2_grp2a3 > apple_noN_w100ks10k.grpAll_grp4_grp3_grp2_grp2a3.avg\n"; 

my $col_cnt = 4; 
my $col_pi  = 5; 
my $col_delt = 3; 

my $h=<>; chomp($h); 
my @ha = split(/\t/, $h); 
my @hh = @ha[0..4]; 
for (my $i=$col_pi; $i<@ha; $i+=$col_delt) {
	push(@hh, "perKb_$ha[$i]"); 
}

print join("\t", @hh)."\n"; 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my @out = @ta[0..4]; 
	for (my $i=$col_pi; $i<@ha; $i+=$col_delt) {
		$ta[$i] //= 0; 
		my $vv = -1; 
		$ta[$col_cnt] > 0 and $vv = $ta[$i] / $ta[$col_cnt] * 1000; 
		push(@out, $vv); 
	}
	print join("\t", @out)."\n"; 
}


