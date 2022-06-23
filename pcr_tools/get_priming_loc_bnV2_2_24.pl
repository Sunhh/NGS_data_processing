#!/usr/bin/perl -w 
# 2010-5-25 11:30:58 input a blastn.o file. 
# 2010-11-03 edit to get clean results. 
# [5/13/2022] Add help information. Add opts
use strict;
use LogInforSunhh;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "min_len:i",  # 15
  "max_diff:i", # 0
  "max_diff_tail:i", # 0
  "max_delt:i", # 6000
  "pri_pen:f",  # 0.1
  "pro_pen:f",  # 1.0
  "help!",
);
$opts{'min_len'}  //= 15;
$opts{'max_diff'} //= 0;
$opts{'max_diff_tail'} //= 0;
$opts{'max_delt'} //= 6000;
$opts{'pri_pen'} //= 0.1;
$opts{'pro_pen'} //= 1.0;

my %primer;
my $min_len  = $opts{'min_len'}; 
my $max_diff = $opts{'max_diff'};
my $max_diff_tail = $opts{'max_diff_tail'}; 
my $max_delt = $opts{'max_delt'}; 
my $pri_pen = $opts{'pri_pen'}; 
my $pro_pen = $opts{'pro_pen'}; 

my %str2num = qw(
+      1
-     -1
1      1
-1    -1
plus   1
minus -1
p      1
m     -1
++     1
+-    -1
); 

#my %perfect_loc; 

my $htxt = <<HH;
####################################################################################################
# perl $0 primer.bn6 > primer.bn6.loc.table
# 
# The primer name format is "(\\S+)_(\\d+)_?[fr]"
# -min_len  =      $min_len
# -max_diff =      $max_diff
# -max_diff_tail = $max_diff_tail
# -max_delt =      $max_delt
# -primer_penalty= $pri_pen
# -pair_penalty  = $pro_pen
HH

!@ARGV and -t and &LogInforSunhh::usage($htxt);



BO: 
while (<>) {
	/^Query name\t|^qseqid\t/ and next; 
	s/\s+$//; 
	my @a = split(/\t/, $_); 
	if ($a[0] =~ /^(\S+)_(\d+)_?([fr])$/i) {
		my ($pro_id, $pri_id, $fr_id) = ($1,$2,$3); 
		# $a[1] == $a[3] or next BO; 
		$a[12] - $a[7] <= $max_diff_tail or next BO; 
		$a[7]-$a[6]+1  >= $min_len or next BO;

		$fr_id = lc($fr_id); 
		#if ($a[3]-$a[2]+1 >= $min_len) {
		# if ($a[1] - ($a[3]-$a[2]+1) <= $max_diff) {
		if ($a[12] - ($a[7]-$a[6]+1) <= $max_diff) {

			# my ($scaf, $scaf_s, $scaf_str) = @a[11, 4, 12]; 
			my ($scaf, $scaf_s, $scaf_str) = @a[1, 8, 14]; 

			defined $str2num{$scaf_str} or die "[Err] Bad strand [$scaf_str]\n"; 
			$scaf_str = $str2num{$scaf_str} * -1; 
			$primer{$pro_id}[$pri_id]{$fr_id}{pen} += $pri_pen; 
			defined $primer{$pro_id}[$pri_id]{$fr_id}{loc}{$scaf} or push(@{$primer{$pro_id}[$pri_id]{$fr_id}{array}}, $scaf); 
			push(@{$primer{$pro_id}[$pri_id]{$fr_id}{loc}{$scaf}}, [$scaf_str, $scaf_s]); 
			#$a[1] == $a[3] and $a[10] == 100 and push(@{$perfect_loc{$pro_id}[$pri_id]{$fr_id}}, [$scaf, $scaf_str, $scaf_s]); 
		}
	}else{
		die "$a[0]\n"; 
	}
}

#print STDERR "file ok.\n"; 

# get perfect 
#

# print STDERR join("\t", 'PrimerPair', 'L_penalty', 'R_penalty', 'P_penalty', 'Est.')."\n"; 

print STDOUT join("\t", qw/PrimerPair targetID targetS targetE targetSize LPenalty RPenalty PPenalty EstPenalty/)."\n";

for my $k1 (sort keys %primer) {
	for (my $i=0; $i<@{$primer{$k1}}; $i++) {
		defined $primer{$k1}[$i] or next; 
		if (defined $primer{$k1}[$i]{f} and defined $primer{$k1}[$i]{r}) { 
			for my $k3 (@{$primer{$k1}[$i]{f}{array}}) { 
				my @ta = @{$primer{$k1}[$i]{f}{loc}{$k3}}; 
				defined $primer{$k1}[$i]{r}{loc}{$k3} or next; 
				my @tb = @{$primer{$k1}[$i]{r}{loc}{$k3}}; 
				for (my $j=0; $j<@ta; $j++) {
					for (my $m=0; $m<@tb; $m++) {
						if ($ta[$j][0] != $tb[$m][0]) {
							my $tt = $ta[$j][0]*$ta[$j][1]+$tb[$m][0]*$tb[$m][1]+1; 
							if ($tt >= 0 and $tt <= $max_delt) {
								$primer{$k1}[$i]{pair}{pen} += $pro_pen; 
								my @aa = sort {$a<=>$b} ($ta[$j][1], $tb[$m][1]); 
								push(@{$primer{$k1}[$i]{pair}{loc}}, [$k3, @aa, $tt]); 
							}
						}
					}#end for m
				}#end for j 
			}#end for k3 
		}
		if (defined $primer{$k1}[$i]{f}) { 
			for my $k3 (@{$primer{$k1}[$i]{f}{array}}) { 
				my @ta = @{$primer{$k1}[$i]{f}{loc}{$k3}}; 
				for (my $j=0; $j<@ta; $j++) {
					for (my $m=$j+1; $m<@ta; $m++) { 
						$m == $j and next; 
						if ($ta[$j][0] != $ta[$m][0]) { 
							my $tt = $ta[$j][0]*$ta[$j][1]+$ta[$m][0]*$ta[$m][1]+1; 
							if ($tt >= 0 and $tt <= $max_delt) {
								$primer{$k1}[$i]{pair}{pen} += $pro_pen; 
								my @aa = sort {$a<=>$b} ($ta[$j][1], $ta[$m][1]); 
								push(@{$primer{$k1}[$i]{pair}{loc}}, [$k3, @aa, $tt]); 
							}
						}
					}#end for m
				}#end for j
			}
		}
		if (defined $primer{$k1}[$i]{r}) {
			for my $k3 (@{$primer{$k1}[$i]{r}{array}}) { 
				my @tb = @{$primer{$k1}[$i]{r}{loc}{$k3}}; 
				for (my $j=0; $j<@tb; $j++) {
					for (my $m=$j+1; $m<@tb; $m++) { 
						$j == $m and next; 
						if ($tb[$j][0] != $tb[$m][0]) {
							my $tt = $tb[$j][0]*$tb[$j][1]+$tb[$m][0]*$tb[$m][1]+1; 
							if ($tt >= 0 and $tt <= $max_delt) {
								$primer{$k1}[$i]{pair}{pen} += $pro_pen; 
								my @aa = sort {$a<=>$b} ($tb[$j][1], $tb[$m][1]); 
								push(@{$primer{$k1}[$i]{pair}{loc}}, [$k3, @aa, $tt]); 
							}
						}
					}#end for m
				}#end for j 
			}
		}
		$primer{$k1}[$i]{f}{pen} -= $pri_pen; 
		$primer{$k1}[$i]{r}{pen} -= $pri_pen; 
		$primer{$k1}[$i]{pair}{pen} -= $pro_pen; 
		my $t1 = $primer{$k1}[$i]{f}{pen}; 
		my $t2 = $primer{$k1}[$i]{r}{pen}; 
		my $t3 = $primer{$k1}[$i]{pair}{pen}; 
		
		if (defined $primer{$k1}[$i]{pair}{loc}) {
			for my $ttr (@{$primer{$k1}[$i]{pair}{loc}}) {
				print join("\t", "${k1}_${i}", @$ttr, $t1, $t2, $t3, sqrt($t1**2+$t2**2)+$t3*10)."\n"; 
			}
		}else{
			print STDERR "No paired location for primer_id [${k1}_$i] found!\n"; 
		}
		
#		print STDERR join ("\t", "${k1}_${i}", $t1, $t2, $t3, sqrt($t1**2+$t2**2)+$t3*10)."\n"; 
#		if ($primer{$k1}[$i]{pair}{pen} >= 1 or ($primer{$k1}[$i]{f}{pen}+$primer{$k1}[$i]{r}{pen})>=1) { 
#			print join("\t", "${k1}_${i}f", $primer{$k1}[$i]{f}{pen}, $primer{$k1}[$i]{pair}{pen}) . "\n"; 
#			print join("\t", "${k1}_${i}r", $primer{$k1}[$i]{r}{pen}, $primer{$k1}[$i]{pair}{pen}) . "\n"; 
#			if (defined $primer{$k1}[$i]{pair}{loc}) {
#				for my $ttr (@{$primer{$k1}[$i]{pair}{loc}}) {
#					print join("\t","", @$ttr)."\n"; 
#				}
#			}
#		}
	}#end for i 
}#end for k1



