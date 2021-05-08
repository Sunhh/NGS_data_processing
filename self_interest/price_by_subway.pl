#!/usr/bin/perl
# Count the cost I spent every month by taking subway for work. 
use strict; 
use warnings; 

my ($cost0, $cost1) = (0,0); 
for (my $d=1; $d<=22; $d++) {
	my ($tc0, $tc1) = &price(7, $cost0, 0.2, 1); 
	$cost0 += $tc0; 
	$cost1 += $tc1; 
	my $diff = sprintf("%0.2f", $cost0-$cost1); 
	print "d=$d 0: sumCost0=$cost0 sumCost1=$cost1 curCost0=$tc0 curCost1=$tc1 save=$diff\n"; 

	if ($d % 5 == 0) {
		($tc0, $tc1) = &price(6, $cost0, 0.2, 1); 
		$cost0 += $tc0; 
		$cost1 += $tc1; 
		# ($tc0, $tc1) = &price(4, $cost0, 0, 1); 
		# $cost0 += $tc0; 
		# $cost1 += $tc1; 
		$diff = sprintf("%0.2f", $cost0-$cost1); 
		print "d=$d 1 *: sumCost0=$cost0 sumCost1=$cost1 curCost0=$tc0 curCost1=$tc1 save=$diff\n"; 
		next; 
		
	}
	($tc0, $tc1) = &price(7, $cost0, 0.2, 1); 
	$cost0 += $tc0; 
	$cost1 += $tc1; 
	$diff = sprintf("%0.2f", $cost0-$cost1); 
	print "d=$d 1: sumCost0=$cost0 sumCost1=$cost1 curCost0=$tc0 curCost1=$tc1 save=$diff\n"; 
}

sub price {
	my ($baseCost, $prevCumCost, $add_discountR, $add_discountCut) = @_; 
	$add_discountR   = 0.2; 
	$add_discountCut = 1; 
	my ($costWoAdd, $costWiAdd); 
	if      ($prevCumCost < 100) {
		; 
	} elsif ($prevCumCost < 150) {
		$baseCost *= 0.8; 
	} elsif ($prevCumCost < 400) {
		$baseCost *= 0.5; 
	}
	$costWoAdd = $baseCost; 
	my $sub = $costWoAdd*$add_discountR; $sub > $add_discountCut and $sub = $add_discountCut; 
	$costWiAdd = $costWoAdd-$sub; 
	return($costWoAdd, $costWiAdd, $prevCumCost+$costWoAdd); 
}# price() 


