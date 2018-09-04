#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 

-t and !@ARGV and die "perl $0 rdCnt.Both.TPM > rdCnt.Both.TPMmean\n"; 

my $hl_txt = <>; 
chomp($hl_txt); 
my @hl_arr = split(/\t/, $hl_txt); 
my @grps; 
my @outs; 
for (my $i=0; $i<@hl_arr; $i++) {
	$hl_arr[$i] =~ m!^(\S+)_Rep(\d+)$!i or do { &tsmsg("[Wrn] Skip bad ID [$hl_arr[$i]]\n"); push(@outs, $i); next; }; 
	my $grpID = $1; 
	if (@grps == 0) {
		push(@grps, [$grpID, [$i], { $hl_arr[$i]=>1 }]); 
	} else {
		my $is_foundGrp = 0; 
		for my $g0 (@grps) {
			if ($g0->[0] eq $grpID) {
				$is_foundGrp = 1; 
				if (defined $g0->[2]{$hl_arr[$i]}) {
					&tsmsg("[Wrn] Skip repeated sample ID [$hl_arr[$i]]\n"); 
					next; 
				} else {
					push(@{$g0->[1]}, $i); 
					$g0->[2]{$hl_arr[$i]} = 1; 
				}
			}
			$is_foundGrp == 1 and last; 
		}
		if ($is_foundGrp == 0) {
			push(@grps, [$grpID, [$i], { $hl_arr[$i]=>1 }]); 
		}
	}
}
@grps > 0 or &stopErr("[Err] No good sample ID found!\n"); 
my @outHeader = @hl_arr[@outs]; 
for my $g0 (@grps) {
	my $grpID = "$g0->[0]_x" . scalar(@{$g0->[1]}); 
	push(@outHeader, $grpID); 
}
print STDOUT join("\t", @outHeader)."\n"; 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my @outLine = @ta[@outs]; 
	my ($sum, $cnt)=(0,0); 
	for my $g0 (@grps) {
		for my $i (@{$g0->[1]}) {
			$sum += $ta[$i]; 
			$cnt ++; 
		}
		push(@outLine, $sum/$cnt); 
	}
	print STDOUT join("\t", @outLine)."\n"; 
}

