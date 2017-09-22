#!/usr/bin/perl 
use strict; 
use warnings; 

use LogInforSunhh; 

-t and !@ARGV and die "perl $0 in.SNP.effect > in.SNP.effect.modest\n"; 

my %order = qw(
  Failed   0
  1.3      1
  1.1      2
  1.2      3 
  2.2      4
  2.1      5
  3.8      6
  3.1      7
  3.2      7
  3.3      7
  3.4      7
  3.5      7
  3.6      7
  3.7      7
); 

my %count; 
while (<>) {
	$. % 100e3 == 1 and &tsmsg("[Msg] $. lines.\n"); 
	chomp; 
	m/^SNP_Effe/i and do { print STDOUT "SNP_type\tSNP_EffMod\t$_\n"; next; }; 
	m/^(\S+)/ or next; 
	my $ee = $1; 
	my @ta1 = split(/,/, $ee); 
	my @class; 
	my $is_failed = 0; 
	my $special_type = ''; 
	if ($ta1[0] !~ m/:/) {
		$special_type = shift(@ta1); 
	}
	for my $tb (@ta1) {
		my @tc = split(/:/, $tb); 
		defined $tc[1] or die "tb=|$tb|\n"; 
		$tc[1] eq 'Failed' and $is_failed = 1; 
		defined $order{$tc[1]} or &stopErr("[Err] Unknown type [$tc[1]]\n"); 
		push(@class, [$tc[1], $tb]); 
	}
	@class = sort { $order{$b->[0]} <=> $order{$a->[0]} } @class; 
	my $modest; 
	my $type; 
	if ( $is_failed == 1 ) {
		if ( $class[0][0] =~ m/^3\./ ) {
			$count{'Failed_Coding'} ++; 
			$type='Failed_Coding'; 
			$modest = "$class[0][1]:Failed_Coding"; 
		} else {
			$count{$class[0][0]} ++; 
			$modest = $class[0][1]; 
			$type = $class[0][0]; 
		}
	} elsif ( $special_type ne '' ) {
		if ( $class[0][0] =~ m/^3\./ ) {
			$type = "${special_type}_coding"; 
			$count{$type} ++; 
			$modest = "$class[0][1]:$type"; 
		} else {
			$type = $class[0][0]; 
			$count{$type} ++; 
			$modest = $class[0][1]; 
		}
	} else {
		$count{$class[0][0]} ++; 
		$modest = $class[0][1]; 
		$type = $class[0][0]; 
	}
	print STDOUT "$type\t$modest\t$_\n"; 
}

=cut 
print STDOUT join("\t", qw/Type Count/)."\n"; 
my ($sumV, $sumC) = (0,0); 
for my $t (sort keys %count) {
	print STDOUT join("\t", $t, $count{$t})."\n"; 
	$sumV += $count{$t}; 
	$sumC ++; 
}
print STDOUT join("\t", "SUM", $sumV)."\n"; 



