#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 
use Parallel::ForkManager; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"refChr:s", # refChr 
	  "refChr_split!", 
	"1col_list:s", # file name of .1col file ; 
	"within_dist:i", # 5 bp
	"cpuN:i", # 1 
	"help!", 
); 

$opts{'within_dist'} //= 5; 
$opts{'cpuN'} //= 1; 

my $help_txt = <<HH; 
################################################################################
# perl $0 -refChr refChr   -1col_list list_1col_filenames   -within_dist 5
#
# -help
#
#   Output .maskClose files. 
#
# -refChr_split    [Boolean]
################################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
for my $tk (qw/refChr 1col_list/) {
	defined $opts{$tk} or &LogInforSunhh::usage($help_txt); 
}
my %IUPAC_b2d; 
{
my @ta = qw(
AA   CC   GG   TT   UU
WAT  SCG  MAC  KGT  RAG   YCT
BCGT DAGT HACT VACG -ACGT NACGT
);
	for my $tb (@ta) {
		my @bb = split(//, $tb);
		my $dbase = shift(@bb);
		$IUPAC_b2d{$dbase} = $dbase;
		for (&permutations(\@bb, scalar(@bb))) {
			$IUPAC_b2d{join('', @$_)} = $dbase;
		}
	}
}# End of IUPAC_xxx;




&tsmsg("[Rec] Begin.\n"); 
my @lis_fn_1col = &_get_1col_list($opts{'1col_list'}); 

if ( ! $opts{'refChr_split'} ) {
	my $refChr = &load_refChr($opts{'refChr'}) ; 
	my $MAX_PROCESSES = $opts{'cpuN'} ; 
	my $pm = new Parallel::ForkManager($MAX_PROCESSES); 
	for my $fn_1col (@lis_fn_1col) {
		my $pid = $pm->start and next; 
		&mask_1colSNP( 'within_dist' => $opts{'within_dist'} , 'fn_1col' => $fn_1col, 'fn_out' => "$fn_1col.maskClose", 'refChr_href'=>$refChr ); 
		$pm->finish; 
	}
	$pm->wait_all_children; 

} else {
	my $tmp_dir = &fileSunhh::new_tmp_dir(); 
	mkdir($tmp_dir) or &stopErr("[Err] mk $tmp_dir/\n"); 
	my %splitted_refChr = &split_refChr( $opts{'refChr'}, "$tmp_dir/refChr" ); 
	my $MAX_PROCESSES = $opts{'cpuN'} ; 
	my $pm = new Parallel::ForkManager($MAX_PROCESSES); 
	for (my $i=0; $i<@lis_fn_1col; $i++) {
		my $pid = $pm->start and next; 
		&split_1col( $lis_fn_1col[$i], "$tmp_dir/$i", \%splitted_refChr ); 
		$pm->finish; 
	}
	$pm->wait_all_children; 
	my @sep_fn_1col_msk; 
	for ( my $i=0; $i<@{$splitted_refChr{'fn'}}; $i++ ) {
		my $refChr = &load_refChr( $splitted_refChr{'fn'}[$i] ); 
		$pm = new Parallel::ForkManager($MAX_PROCESSES); 
		for ( my $j=0; $j<@lis_fn_1col; $j++ ) {
			my $fn_1col = "$tmp_dir/$j.$i"; 
			push(@{$sep_fn_1col_msk[$j]}, "$fn_1col.maskClose"); 
			my $pid = $pm->start and next; 
			&mask_1colSNP( 'within_dist' => $opts{'within_dist'} , 'fn_1col' => $fn_1col, 'fn_out' => "$fn_1col.maskClose", 'refChr_href'=>$refChr ); 
			$pm->finish; 
		}
		$pm->wait_all_children; 
	}
	for (my $i=0; $i<@lis_fn_1col; $i++) {
		&cat_1col("$lis_fn_1col[$i].maskClose", $sep_fn_1col_msk[$i] ); 
	}
	&fileSunhh::_rmtree($tmp_dir); 
}
&tsmsg("[Rec] All done.\n"); 

sub cat_1col {
	my ($out_fn, $in_fn_aref) = @_; 
	my $o_fh = &openFH($out_fn, '>'); 
	for (my $i=0; $i<@$in_fn_aref; $i++) {
		my $in_fh = &openFH($in_fn_aref->[$i], '<'); 
		my $hh = readline($in_fh); 
		$i == 0 and print {$o_fh} $hh; 
		while (<$in_fh>) {
			print {$o_fh} $_; 
		}
	}
	return; 
}# cat_1col () 

sub split_1col {
	my ($fn, $opref, $spl) = @_; 
	my $fh = &openFH($fn,'<'); 
	my $hh = readline($fh); 
	my $fh_o; 
	while (<$fh>) {
		if (defined $spl->{'ln2idx'}{$.}) {
			defined $fh_o and close($fh_o); 
			my $o_fn = "${opref}.$spl->{'ln2idx'}{$.}"; 
			$fh_o = &openFH($o_fn, '>'); 
			print {$fh_o} $hh; 
		}
		print {$fh_o} $_; 
	}
	return; 
}# split_1col 

sub split_refChr {
	my ($fn, $opref) = @_; 
	my %back; 
	my $fh = &openFH($fn,'<'); 
	my $fh_o; 
	my $hh = readline($fh); 
	my $prevID = ''; 
	$back{'idx'} = -1; 
	while (<$fh>) {
		$_ =~ m!^(\S+)\t(\d+)\t(\S+)$! or die "bad line: $_\n"; 
		if ( $prevID ne $1 ) {
			$prevID = $1; 
			&tsmsg("[Msg] Generating splitted refChr [$fn] [$prevID]\n"); 
			$back{'idx'} ++; 
			my $o_fn = "${opref}.$back{'idx'}"; 
			$back{'ln2fn'}{$.} = $o_fn; 
			$back{'ln2idx'}{$.} = $back{'idx'}; 
			push(@{$back{'fn'}}, $o_fn); 
			defined $fh_o and close($fh_o); 
			$fh_o = &openFH($o_fn, '>'); 
			print {$fh_o} $hh; 
		}
		print {$fh_o} $_; 
	}
	close($fh); 
	close($fh_o); 
	return(%back); 
}

# Required : 'fn_1col/fh_i' , 'refChr_href'; 
sub mask_1colSNP {
	my %parm = &mathSunhh::_setHashFromArr( @_ ); 
	$parm{'within_dist'}   //= 5; 
	$parm{'fn_out'}        //= "$parm{'fn_1col'}.maskClose"; 
	$parm{'fn_outType'}    //= '>'; 
	$parm{'close_fh_i'}    //= ( (defined $parm{'fh_i'}) ? 0 : 1 ); 
	$parm{'close_fh_o'}    //= ( (defined $parm{'fh_o'}) ? 0 : 1 ); 
	$parm{'fh_i'}          //= &openFH($parm{'fn_1col'}, '<'); 
	$parm{'fh_o'}          //= &openFH($parm{'fn_out'}, $parm{'fn_outType'}); 
	$parm{'need_header'}   //= 1; 
	$parm{'out_tail'}      //= 1; 
	if ( $parm{'need_header'} ) {
		my $head_line = readline($parm{'fh_i'}); 
		print {$parm{'fh_o'}} $head_line; 
	}

	my %wind; # Keys : qw/siteN wind_start_idx curr_base_idx wind_bases idx2foll/
	if ( defined $parm{'wind'} ) {
		%wind = %{$parm{'wind'}}; 
	} else {
		&reset_wind( \%wind ); 
		$wind{'curr_base_idx'} = -1; 
	}
	my @windIdx_to_follow; # $windIdx_to_follow[0] = [0..4]; $windIdx_to_follow[1] = [1,2,3,4,0]; ... 
	for ( my $i=0; $i<$parm{'within_dist'}; $i++ ) {
		for ( my $j=$i; ; $j++) {
			$j >= $parm{'within_dist'} and $j -= $parm{'within_dist'}; 
			if ( $i == 0 ) {
				push( @{$windIdx_to_follow[$i]}, $j ); 
				$j == $parm{'within_dist'} - 1 and last; 
			} else {
				push( @{$windIdx_to_follow[$i]}, $j ); 
				$j == $i-1 and last; 
			}
		}
	}
	$wind{'idx2foll'} = \@windIdx_to_follow; 

	my %tmp_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>5e6 ); 
	for my $ref_base ( @{ $parm{'refChr_href'}{'base'} } ) {
		$_ = readline($parm{'fh_i'}); 
		chomp; 
		&fileSunhh::log_section($. , \%tmp_cnt) and &tsmsg("[Msg] Processing $. line in [$parm{'fn_1col'}]. $_\n"); 
		$wind{'curr_base_idx'} ++; 
		# $. == $wind{'curr_base_idx'}+2 or &stopErr("[Err] line number not match. $. == $wind{'curr_base_idx'}+2\n"); 

		my $base_compare = &compare_base_1col( $_, $ref_base ); 

		if ( defined $parm{'refChr_href'}{'split_idx'}{ $wind{'curr_base_idx'} } and $wind{'siteN'} > 0 ) {
			# Clear and reset %wind. 
			my @all_bp = &all_first_bp( \%wind ); 
			for my $first_bp (@all_bp) {
				print {$parm{'fh_o'}} "$first_bp\n"; 
			}
			&reset_wind( \%wind ); 
		} elsif ( $wind{'siteN'} == $parm{'within_dist'} ) {
			my $first_bp = &new_first_bp( \%wind ); 
			print {$parm{'fh_o'}} "$first_bp\n"; 
		} elsif ( $wind{'siteN'} < $parm{'within_dist'} ) {
		} else {
			&stopErr("[Err] bad siteN [$wind{'siteN'}]\n"); 
		}
		&update_wind( \%wind, $base_compare, $parm{'within_dist'} ); 

	}
	if ( $parm{'out_tail'} ) {
		for my $first_bp ( &all_first_bp( \%wind ) ) {
			print {$parm{'fh_o'}} "$first_bp\n"; 
		}
		&reset_wind(\%wind); 
	}
	$parm{'close_fh_o'} and close($parm{'fh_o'}); 
	$parm{'close_fh_i'} and close($parm{'fh_i'}); 
	return( \%wind ); 
}#mask_1colSNP () 

sub reset_wind {
	$_[0]->{'siteN'} = 0; 
	$_[0]->{'wind_start_idx'} = 0; 
	$_[0]->{'wind_bases'} = []; 
}# 
sub all_first_bp {
	# ($wind_href) = @_; 
	my $curr_idx = $_[0]->{'wind_start_idx'}; 
	$_[0]->{'siteN'} > 0 or return; 
	my $stop_idx = $curr_idx - 1; $stop_idx == -1 and $stop_idx = $_[0]->{'siteN'} - 1; 
	my @back; 
	for (my $i=$curr_idx; ; $i++) {
		$i >= $_[0]->{'siteN'} and $i -= $_[0]->{'siteN'}; 
		$_[0]->{'wind_start_idx'} = $i; 
		my $first_bp = &new_first_bp( $_[0] ); 
		push(@back, $first_bp); 
		$i == $stop_idx and last; 
	}
	return(@back); 
}# all_first_bp () 

# Updated : qw/siteN wind_bases wind_start_idx/ 
# $wind{'wind_bases'}[$i] = [ $base_1col, $diff2ref(1,2,3,4), $mask(1,0) ]
sub update_wind {
	# ( $wind_href , $base_compare, $within_dist ) = ( \%wind, $base_compare, $parm{'within_dist'} ); 
	if ( $_[0]->{'siteN'} == $_[2] ) {
		# Rolling. 
		# 'wind_start_idx' is to be changed. 

		my $mask = 0; 
		my $curr_idx = $_[0]->{'wind_start_idx'}; 
		my @foll = @{ $_[0]->{'idx2foll'}[$curr_idx] }; 
		if ( $_[1][1] == 1 ) {
			for ( my $i=1; $i<@foll; $i++ ) {
				if ( $_[0]->{'wind_bases'}[ $foll[$i] ][1] == 1 ) {
					$mask = 1; 
					for ( my $j=$i; $j<@foll; $j++ ) {
						$_[0]->{'wind_bases'}[ $foll[$j] ][2] = 1; 
					}
					last; 
				}
			}
		}
		$_[0]->{'wind_bases'}[ $curr_idx ] = [ $_[1][0], $_[1][1], $mask ]; 
		$_[0]->{'wind_start_idx'} = $foll[1]; 

	} elsif ( $_[0]->{'siteN'} < $_[2] ) {
		$_[0]->{'wind_start_idx'} = 0; 
		$_[0]->{'siteN'} ++; 
		my $mask = 0; 
		if ( $_[1][1] == 1 ) {
			# Add a different sites; 
			for ( my $i=0; $i<@{ $_[0]->{'wind_bases'} }; $i++ ) {
				if ( $_[0]->{'wind_bases'}[$i][1] == 1 ) {
					$mask = 1; 
					for ( my $j=$i; $j<@{ $_[0]->{'wind_bases'} }; $j++ ) {
						$_[0]->{'wind_bases'}[$i][2] = 1; 
					}
					last; 
				}
			}
		}
		push( @{ $_[0]->{'wind_bases'} }, [ $_[1][0], $_[1][1], $mask ] ); 
	} else {
		&stopErr("[Err] Too big siteN [$_[0]->{'siteN'}]\n"); 
	}
	return; 
}# update_wind() 

=head1 compare_base_1col( $base_1col, uc($base_ref) ) 

Return : ( [ $base_1col, $diff:1 or $same:2 or $ref_as_N:3 or $1col_as_N:4 ] )

=cut
sub compare_base_1col {
	# Input : ( $base_1col , $base_ref ) = @_; 
	my $base_1col = &cnvt_base( $_[0] ); 
	if ( $base_1col eq 'N' ) {
		return( [ $_[0], 4 ] ); 
	} elsif ( $_[1] eq 'N' ) {
		return( [ $_[0], 3 ] ); 
	} elsif ( $base_1col eq $_[1] ) {
		return( [ $_[0], 2 ] ); 
	} else {
		return( [ $_[0], 1 ] ); 
	}
	&stopErr("[Err] Wrong get here.\n"); 
}# compare_base_1col() 

sub cnvt_base {
	# I am thinking to make it as degenerated base. 
	my $tb = uc($_[0]); 
	defined $IUPAC_b2d{$tb} and $tb = $IUPAC_b2d{$tb}; 
	return( $tb ); 
}

=head1 load_refChr( 'refChr' )

Return       : ( \%refChr )

 $refChr{'base'} = [ $base_pos1, $base_pos2, ... ]; # Skipping the first line. All bases are upper cases. 
 $refChr{'split_idx'}{ $split_base_pos } = $line_number; 
   At $split_base_pos, the previous position has different chrID from the current chrID. 

File format of refChr : 

#   chr	pos	base
#   chr10	1	T
#   chr10	2	G
#   chr10	3	T

=cut
sub load_refChr {
	my $fn = shift @_; 
	my $fh = &openFH($fn,'<'); 
	my %back; 
	$back{'base_idx'} = -1; 
	my %tmp_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>5e6 ); 
	<$fh>; 
	while (<$fh>) {
		&fileSunhh::log_section( $. , \%tmp_cnt ) and &tsmsg("[Msg] Loading $. line [$fn].\n"); 
		chomp; 
		$_ =~ m!^(\S+)\t(\d+)\t(\S+)$! or die "bad line: $_\n"; 
		# my ($cID, $pos, $base) = ($1, $2, $3); 
		$back{'prevChr'} //= $1; 
		push(@{$back{'base'}}, uc($3)); 
		$back{'base_idx'} ++; 
		if ( $back{'prevChr'} ne $1 ) {
			$back{'prevChr'} = $1; 
			$back{'split_idx'}{ $back{'base_idx'} } = $.; 
			&tsmsg("[Msg] Switching chrID at line : $_\n"); 
		}
	}
	close($fh); 
	delete($back{'prevChr'}); 
	&tsmsg("[Msg] refChr [$fn] loaded.\n"); 
	return(\%back); 
}# load_refChr() 

sub _get_1col_list {
	my $fn = shift(@_); 
	my $fh = &openFH($fn, '<'); 
	my @back; 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		push(@back, $ta[0]); 
	}
	close($fh); 
	return(@back); 
}# _get_1col_list () 

sub new_first_bp {
	# ( $wind_href ) = ( \%wind ); 
	if ( $_[0]->{'wind_bases'}[ $_[0]->{'wind_start_idx'} ][ 2 ] == 1 ) {
		return('N'); 
	} else {
		return( $_[0]->{'wind_bases'}[ $_[0]->{'wind_start_idx'} ][ 0 ] ); 
	}
}# new_first_bp () 

sub permutations {
	my ($list, $n) = @_;
	$n = $n // scalar(@$list);
	$n > @$list and return ($list);
	$n <= 1 and return(map {[$_]} @$list);
	my @perm;
	for my $i (0 .. $#$list) {
		my @rest = @$list;
		my $val = splice(@rest, $i, 1);
		for ( &permutations(\@rest, $n-1) ) {
			push(@perm, [$val, @$_]);
		}
	}
	return @perm;
}#sub permutations()

