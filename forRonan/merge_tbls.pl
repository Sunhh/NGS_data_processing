#!/usr/bin/perl
use strict; 
use warnings; 
use mathSunhh; 
use fileSunhh; 
use LogInforSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"basicTbl:s", # Format : [ID, Start, End, ...]
	"addTbl:s@", # Format : [ID, Start, End]
	"addTbl_list:s", 
	"detail!", 
	"noValue!", 
	"outfile:s", 
); 

sub usage {
	print <<HH; 
################################################################################
# perl $0 -basicTbl Spec97_loc.tbl -addTbl A1_pileup.tbl -addTbl A2_pileup.tbl
# 
# -help 
# -detail 
# 
# -basicTbl           The basic table includes all positions that will be considered. 
# -addTbl             The position tables that will be attached to basic table. 
# -addTbl_list        This will overwrite -addTbl parameter. 
# 
# -noValue            [Boolean] If given, No value adding step processed. 
# -outfile            [filename] Output to STDOUT if not given. 
# 
################################################################################
HH
	exit 1; 
}

$opts{'help'} and &usage(); 
defined $opts{'basicTbl'} or &usage(); 
defined $opts{'addTbl'} or defined $opts{'addTbl_list'} or &usage(); 
if (defined $opts{'addTbl_list'}) {
	my $alFh = &openFH($opts{'addTbl_list'}, '<'); 
	@{$opts{'addTbl'}} = (); 
	while (<$alFh>) {
		chomp; 
		m/^\s*(#|$)/ and next; 
		my @ta = split(/\s/, $_); 
		push(@{$opts{'addTbl'}}, $ta[0]); 
	}
}
$opts{'outFh'} = \*STDOUT; 
defined $opts{'outfile'} and $opts{'outFh'} = &openFH( $opts{'outfile'}, '>' ); 

my %basic_tbl; 
my @basic_IDs; 

&tsmsg("[Rec] Beginning $0.\n"); 

{
	&tsmsg("[Msg] Reading basic table [$opts{'basicTbl'}]\n"); 
	my $btFh = &openFH($opts{'basicTbl'}, '<'); 
	my ($t_h, $t_a) = &tblFh2sorted_lociHash( $btFh, 1 ); 
	%basic_tbl = %$t_h; 
	@basic_IDs = @$t_a; 
	close($btFh); 
}
$opts{'detail'} and &tsmsg("[Msg] Total ", scalar(@basic_IDs) ," scaffolds sorted.\n"); 

&tsmsg("[Msg] Generating table bins.\n"); 
for my $at (@{$opts{'addTbl'}}) {
	&tsmsg("[Msg] Adding [$at]\n"); 
	my $atFh = &openFH($at, '<'); 
	my ($t_h, $t_a) = &tblFh2sorted_lociHash( $atFh, 1 ); 
	close($atFh); 
	
	my %ovlpSE_add; 
	for my $add_ID (@$t_a) {
		defined $basic_tbl{$add_ID} or next; 
		my $idx_start = 0; 
		my $basic_tbl_maxID = @{$basic_tbl{$add_ID}} - 1; 
		for my $addSE (@{$t_h->{$add_ID}}) {
			my ( $specSE, $ovlpSE, $idx_start_new ) = &mathSunhh::sep_loc2_by_loc1_singleLoc2( 
			  [ @{$basic_tbl{$add_ID}}[ $idx_start .. $basic_tbl_maxID ] ], 
			  [ $addSE->[0], $addSE->[1] ]
			); 
			$idx_start_new > 0 and $idx_start_new -- ; # When the add_loci is included by the previous basic_loci. 
			$idx_start += $idx_start_new ;             # Convert the index back to $basic_tbl{$add_ID} ; 
			$idx_start > $basic_tbl_maxID and &stopErr("[Err] idx_start=$idx_start, basic_tbl_maxID=$basic_tbl_maxID, add_ID=$add_ID");    # This should never happen. 
			@$ovlpSE > 0 or next;                      # No overlapping region betwee added and basic. 
			push(@{$ovlpSE_add{$add_ID}}, @$ovlpSE); 
		}
	}
	for my $add_ID (keys %ovlpSE_add) {
		@{$ovlpSE_add{$add_ID}} = sort { $a->[0]<=>$b->[0] || $a->[1]<=>$b->[1] } @{$ovlpSE_add{$add_ID}}; 
	}
	
	&tsmsg("[Msg]   Splitting tables.\n"); 
	my %backSE_new; 
	for my $id (@basic_IDs) {
		unless ( defined $ovlpSE_add{$id} ) {
			$backSE_new{$id} = [ @{$basic_tbl{$id}} ]; 
			next; 
		}
		my $idx_start = 0; 
		my $ovlpSE_add_maxID = @{$ovlpSE_add{$id}} - 1; 
		for my $basicSE (@{$basic_tbl{$id}}) {
			my ( $specSE, $ovlpSE,  $idx_start_new ) = &mathSunhh::sep_loc2_by_loc1_singleLoc2( 
			  [ @{$ovlpSE_add{$id}}[ $idx_start .. $ovlpSE_add_maxID ] ], 
			  [$basicSE->[0], $basicSE->[1]] 
			); 
			$idx_start_new > 0 and $idx_start_new --; 
			$idx_start += $idx_start_new; 
			$idx_start > $ovlpSE_add_maxID and &stopErr("[Err] idx_start=$idx_start, maxID=$ovlpSE_add_maxID, basicID=$id\n"); 
			my @a = sort { $a->[0]<=>$b->[0] || $a->[1]<=>$b->[1] } ( @$specSE, @$ovlpSE); 
			for my $ar1 (@a) {
				push( @{$backSE_new{$id}}, [ $ar1->[0], $ar1->[1] ] ); 
			}
		}
	}
	%basic_tbl = %backSE_new; 
}

unless ($opts{'noValue'}) {
	&tsmsg("[Msg] Adding table values in each bin.\n"); 
	my $idx_add = -1; 
	for my $at ( @{$opts{'addTbl'}} ) {
		$idx_add ++; 
		&tsmsg("[Msg] Adding $idx_add file [$at]\n"); 
		my $atFh = &openFH($at, '<'); 
		while (<$atFh>) {
			$opts{'detail'} and $. % 100e3 == 1 and &tsmsg("[Msg]   Processing $. line\n"); 
			chomp; 
			my @ta = split(/\t/, $_); 
			my ($id, $s, $e) = @ta; 
			defined $basic_tbl{$id} or next; 
			my ( $specSE, $ovlpSE ) = &mathSunhh::sep_loc2_by_loc1_singleLoc2( $basic_tbl{$id}, [$s, $e] ); 
			@$ovlpSE > 0 or next; 
			for my $ar1 (@$ovlpSE) {
			#&tsmsg("ar1: @$ar1\n"); 
				for my $ar2 ( @{$basic_tbl{$id}} ) {
				#&tsmsg("ar2: @$ar2\n"); 
					if ( $ar2->[0] == $ar1->[0] and $ar2->[1] == $ar1->[1] ) {
						$ar2->[$idx_add + 2] = $ar1->[1]-$ar1->[0]+1; 
					} else {
						$ar2->[$idx_add + 2] = -1; 
					}
				}
			}
		}
		close($atFh); 
		for my $id (@basic_IDs) {
			for my $ar2 (@{$basic_tbl{$id}}) {
				$ar2->[$idx_add + 2] //= -1; 
			}
		}
	}
}#End unless ($opts{'noValue'}) 

&tsmsg("[Rec] Output result.\n"); 
if ($opts{'noValue'}) {
	print {$opts{'outFh'}} join("\t", qw/RefID RefStart RefEnd/)."\n"; 
} else {
	print {$opts{'outFh'}} join("\t", qw/RefID RefStart RefEnd/, @{$opts{'addTbl'}})."\n"; 
}

for my $id ( @basic_IDs ) {
	for my $ar1 (@{$basic_tbl{$id}}) {
		print {$opts{'outFh'}} join("\t", $id, @$ar1)."\n"; 
	}
}

&tsmsg("[Rec] All done for $0\n"); 

####################### Sub-routines. 
sub tblFh2sorted_lociHash{
	my ($tblFh, $ifSort) = @_;  
	$ifSort //= 1; 
	my (%back_hash, @back_ID); 
	while (<$tblFh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		my ($id, $s, $e) = @ta[0,1,2]; 
		$id eq 'RefID' and next; 
		defined $back_hash{$id} or push(@back_ID, $id); 
		push(@{$back_hash{$id}}, [$s, $e]); 
	}
	if ($ifSort) {
		for my $id (@back_ID) {
			@{$back_hash{$id}} = sort { $a->[0]<=>$b->[0] || $a->[1]<=>$b->[1] } @{$back_hash{$id}}; 
		}
	}
	return (\%back_hash, \@back_ID); 
}# tblFh2sorted_lociHash() 


