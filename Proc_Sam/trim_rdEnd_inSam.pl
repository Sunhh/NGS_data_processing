#!/usr/bin/perl
# Please note that the mate-mapping information will be wrong here. 
use strict; 
use warnings; 
use SeqAlnSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"trimLen:i", # Default 10; 
); 
$opts{'trimLen'} //= 10; 

sub usage {
	print STDERR <<HH; 

perl $0 in.sam > out.sam

-help
-trimLen      [$opts{'trimLen'}] The length to be trimmed from both ends. 

HH
	exit(1); 
}

-t and !@ARGV and &usage(); 
$opts{'help'} and &usage(); 

my %flag_fwd = %{ &SeqAlnSunhh::mk_flag('keep'=>'2=0,4=0') }; 
my %flag_rev = %{ &SeqAlnSunhh::mk_flag('keep'=>'2=0,4=1') }; 


while (<>) {
	chomp; 
	my @ta = (split(/\t/, $_)); 
	if ($ta[0] =~ m/^@/ and @ta <= 10) {
		# Headers. 
		print STDOUT "$_\n"; 
		next; 
	}
	
	unless (defined $flag_fwd{$ta[1]} or defined $flag_rev{$ta[1]}) {
		print STDOUT "$_\n"; 
		next; 
	}

	# rdID  flag  refID  refPos  MapQ  Cigar(100M)   MateRefID MateRefPos Insert ReadBases qualities # Others   AS:i:-15        XN:i:0  XM:i:3  XO:i:0  XG:i:0  NM:i:3  MD:Z:3G43G46A5  YS:i:0  YT:Z:DP
	# 0     1     2      3       4     5             6         7          8      9         10        11
	my ($new_cigar_aref, $new_drop_lp, $new_drop_rp, $left_mv_refP, $right_mv_refP) = &trimCigar($ta[5], $opts{'trimLen'}); 
	scalar( @$new_cigar_aref ) > 0 or next; 
	my $new_cigar_str = join('', map { "$_->[0]$_->[1]" } @$new_cigar_aref); 
	my $rd_len = length($ta[9]); 
	my $new_rd_bp = substr( $ta[9], $new_drop_lp, $rd_len-$new_drop_lp-$new_drop_rp ); 
	my $new_rd_qb = substr( $ta[10], $new_drop_lp, $rd_len-$new_drop_lp-$new_drop_rp ); 
	my $new_refP = $ta[3]+$left_mv_refP; 
	

	print STDOUT join("\t", @ta[0 .. 2], $new_refP, $ta[4], $new_cigar_str, @ta[6 .. 8], $new_rd_bp, $new_rd_qb)."\n"; 
}

sub trimCigar {
	# $_[0] is cigar string. 
	# $_[1] is the end length for trimming. 
	my $left_trimmed_RdP = 0; 
	my $right_trimmed_RdP = 0; 
	my $left_mv_refP = 0; 
	my $right_mv_refP = 0; 
	my @cigar_arr; 

	while ($_[0] =~ s/^(\d+)([MIDNSHP=X])//) {
		my ($len, $tag) = ($1, $2); 
		push(@cigar_arr, [$len, $tag]); 
	}

	# Check from left 
	my @new_cigar_arr; 
	my $trimmed_RdP = 0; 
	my $bp = 0; 
	for (@cigar_arr) {
		if ($trimmed_RdP >= $_[1]) {
			push(@new_cigar_arr, $_); 
			next; 
		}
		if ($_->[1] =~ m/^[MSP=X]$/i) {
			my $raw_rdP = $trimmed_RdP; 
			($trimmed_RdP, $bp) = &trim_by_pos($trimmed_RdP, $_[1], $_->[0]); 
			$bp > 0 and push(@new_cigar_arr, [$bp, $_->[1]]); 
			if ( $_->[1] =~ m/^[MP=X]$/i ) {
				$left_mv_refP += ($trimmed_RdP-$raw_rdP); 
			}
		} elsif ($_->[1] =~ m/^I$/i) {
			$trimmed_RdP += $_->[0]; 
			$bp = 0; 
		} elsif ($_->[1] =~ m/^D$/i) {
			$left_mv_refP += $_->[0]; 
		}
	}
	$left_trimmed_RdP = $trimmed_RdP; 
	# Check from right
	@cigar_arr = reverse(@new_cigar_arr); 
	@new_cigar_arr = (); 
	$trimmed_RdP = $bp = 0; 
	for (@cigar_arr) {
		if ($trimmed_RdP >= $_[1]) {
			push(@new_cigar_arr, $_); 
			next; 
		}
		if ($_->[1] =~ m/^[MSP=X]$/i) {
			my $raw_rdP = $trimmed_RdP; 
			($trimmed_RdP, $bp) = &trim_by_pos($trimmed_RdP, $_[1], $_->[0]); 
			$bp > 0 and push(@new_cigar_arr, [$bp, $_->[1]]); 
			if ( $_->[1] =~ m/^[MP=X]$/i ) {
				$right_mv_refP += ($trimmed_RdP-$raw_rdP); 
			}
		} elsif ($_->[1] =~ m/^I$/i) {
			$trimmed_RdP += $_->[0]; 
			$bp = 0; 
		} elsif ($_->[1] =~ m/^D$/i) {
			$right_mv_refP += $_->[0]; 
		}
	}
	$right_trimmed_RdP = $trimmed_RdP; 

	@cigar_arr = reverse(@new_cigar_arr); 

	return (\@cigar_arr, $left_trimmed_RdP, $right_trimmed_RdP, $left_mv_refP, $right_mv_refP); 
	
}

sub trim_by_pos {
	# ( Start_RdP-1, Max_to_trim, len_input)
	if      ( $_[1]-$_[0] < 1 ) {
		# Cannot trimmed more. 
		return ($_[0], $_[2]); 
	} elsif ( $_[1]-$_[0] >= $_[2] ) {
		return ($_[2]+$_[0], 0); 
	} elsif ( $_[1]-$_[0] < $_[2] ) {
		return ( $_[1], $_[2]-($_[1]-$_[0]) ); 
	} else {
		die "why here? [@_]\n"; 
	}
}




