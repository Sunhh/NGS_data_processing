#!/usr/bin/perl -w
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"task:s", # work to do
	"in_blkTxt:s", # blocks_coords.txt ; 
	"in_blkGff:s", # blocks_coords.gff ; 
	"in_scfLen:s", # in.fasta.key_len ; 
	"out_blk:s", # Will be STDOUT if not given. 
	 "scfID_order:s@", # Default is the first scaffold. 
	 "scfID_order_list:s", # A list that defines the order of scfIDs; 
	"chk_ctgStrand:s@", 
); 

$opts{'help'} and &_usage(); 
defined $opts{'task'} or &_usage(); 
$opts{'task'} //= ''; 

my %glob; 

if ($opts{'task'} =~ m!^coords_txt2gff$!i) {
	&cnvt_coords_txt2gff(); 
} elsif ( $opts{'task'} =~ m!^coords_gff2txt$!i ) {
	&cnvt_coords_gff2txt(); 
} elsif ( $opts{'task'} =~ m!^sort_blkGff$!i ) {
	&sort_coords_gff(); 
} elsif ( $opts{'task'} =~ m!^fix_strand_gff$!i ) {
	# If the first segment in defined reference scaffold has '-' strand, reverse all segments in this block; 
	# If there is no reference scaffold in the block, using the first segment in the block. 
	&fix_strand_gff(); 
} elsif ( $opts{'task'} =~ m!^coords_gff2tab$!i ) {
	# Output format : scfID_1 \\t scfStrand_1 \\t scfStart_1 \\t scfEnd_1 \\t scfID_2 \\t scfStrand_2 \\t scfStart_2 \\t scfEnd_2 \\t srcBlkID \\n
	# Only blocks with even number of segments are considered. 
	&cnvt_coords_gff2tab(); 
}


=head1 _usage()

Function    : Output help information and then exit. 

=cut
sub _usage {
	my $help_txt = <<HH; 
################################################################################
# perl $0 -task some_task
#
# -help 
#
# -task 
#   'coords_txt2gff' : Convert Sibelia output txt file to gff format; 
#     Need -in_blkTxt [filename], could be 'blocks_coords.txt'; 
#   'coords_gff2txt' : Convert Sibelia output gff file to txt format; 
#     Need -in_blkGff [filename], could be 'blocks_coords.gff'; 
#   
#   'sort_blkGff'    : Sort input block gff according to scaffolds order and positions. 
#     Need -in_blkGff . 
#     Could use -scfID_order_list , which is a list scaffold IDs, one scfID per line. 
#     Could use " -scfID_order scfID_1 -scfID_order scfID_2 ... " for selected scfIDs. 
#   'fix_strand_gff' : Change the representative segment to '+' strand. 
#     Need -in_blkGff , and could use -scfID_order_list and -scfID_order\@ ; 
#   'coords_gff2tab' : Exchange paired segments if the first segment should be after the 2nd one. 
#     Need -in_blkGff , and could use -scfID_order_list and -scfID_order\@ ;
#     -chk_ctgStrand\@  [agp_filename] If given, two additional columns output to show the strand in segments. 
################################################################################
HH
	&LogInforSunhh::usage($help_txt); 
	return; 
}# _usage() 

################################################################################
# Using sub-routines. 
################################################################################

sub cnvt_coords_gff2tab {
	# Output format : scfID_1 \\t scfStrand_1 \\t scfStart_1 \\t scfEnd_1 \\t scfID_2 \\t scfStrand_2 \\t scfStart_2 \\t scfEnd_2 \\t srcBlkID \\n
	my $outFh = \*STDOUT; 
	defined $opts{'out_blk'} and $outFh = &openFH( $opts{'out_blk'}, '>' ); 
	my %data = %{ &_load_coords_gff( $opts{'in_blkGff'} ) }; 
	my %new_data; 
	for my $blk_num (keys %{$data{'blk_num2inf'}}) {
		scalar(@{$data{'blk_num2inf'}{$blk_num}}) % 2 == 0 or do { &tsmsg("[Wrn] Skip bad blk_num [$blk_num]\n"); next; }; 
		$new_data{'blk_num2inf'}{$blk_num} = $data{'blk_num2inf'}{$blk_num}; 
		$new_data{'blk_num2order'}{$blk_num} = $data{'blk_num2order'}{$blk_num}; 
	}
	for my $tk (keys %data) {
		defined $new_data{$tk} and next; 
		$new_data{$tk} = $data{$tk}; 
	}
	%data = %new_data; 
	undef(%new_data); 

	&_fix_pairScf_order( \%data ); 
	
	my %scf2ctg; 
	if (defined $opts{'chk_ctgStrand'}) {
		my %has_ctg; 
		for my $fn ( @{$opts{'chk_ctgStrand'}} ) {
			my $c2s = &fileSunhh::load_agpFile($fn, 1); 
			for my $ctgID (keys %$c2s) {
				for my $a1 ( @{$c2s->{$ctgID}} ) {
					push( @{$scf2ctg{$a1->[2]}}, [@{$a1}[3,4], $ctgID, @{$a1}[0,1,5]] ); 
				}
			}
		}
		print {$outFh} join("\t", qw/scf1_ID scf1_Str scf1_Start scf1_End seg1_len scf2_ID scf2_Str scf2_Start scf2_End seg2_len blk_num scf1_ctgStr scf2_ctgStr/)."\n"; 
	} else {
		print {$outFh} join("\t", qw/scf1_ID scf1_Str scf1_Start scf1_End seg1_len scf2_ID scf2_Str scf2_Start scf2_End seg2_len blk_num/)."\n"; 
	}
	for my $blk_num (sort { $data{'blk_num2order'}{$a} <=> $data{'blk_num2order'}{$b} } keys %{$data{'blk_num2inf'}}) {
		my @ta = @{ $data{'blk_num2inf'}{$blk_num} }; 
		for (my $i=0; $i<@ta; $i+=2) {
			my $j=$i+1; 
			if ( defined $opts{'chk_ctgStrand'}) {
				my %cnt; 
				$cnt{'i'}{'scf_info'} = [ @{$ta[$i]}[0,2,3] ]; 
				$cnt{'j'}{'scf_info'} = [ @{$ta[$j]}[0,2,3] ]; 
				
				for my $ij (qw/i j/) {
					for my $t1 (qw/F R U/) {
						$cnt{$ij}{$t1} = 0; 
					}
					for my $a1 ( @{$scf2ctg{$cnt{$ij}{'scf_info'}[0]}} ) {
						$a1->[0] > $cnt{$ij}{'scf_info'}[2] and next; 
						$a1->[1] < $cnt{$ij}{'scf_info'}[1] and next; 
						my ($ovl_len, $ovl_se) = $ms_obj->ovl_region( $cnt{$ij}{'scf_info'}[1],$cnt{$ij}{'scf_info'}[2], $a1->[0], $a1->[1] ); 
						$ovl_len > 0 or &stopErr("[Err] $cnt{$ij}{'scf_info'}[1],$cnt{$ij}{'scf_info'}[2], $a1->[0], $a1->[1]\n"); 
						if ( $a1->[5] eq '+' ) {
							$cnt{$ij}{'F'} += $ovl_len; 
						} elsif ( $a1->[5] eq '-' ) {
							$cnt{$ij}{'R'} += $ovl_len; 
						} elsif ( $a1->[5] eq '?' ) {
							$cnt{$ij}{'U'} += $ovl_len; 
						} else {
							&stopErr("[Err] Unknown strand [$a1->[5]]\n"); 
						}
					}
				}
				my $str_i = join("_", ( $cnt{'i'}{'F'} > 0 or $cnt{'i'}{'R'} > 0 ) ? 'Y' : 'N' , map { "${_}$cnt{'i'}{$_}" } qw/F R U/); 
				my $str_j = join("_", ( $cnt{'j'}{'F'} > 0 or $cnt{'j'}{'R'} > 0 ) ? 'Y' : 'N', map { "${_}$cnt{'j'}{$_}" } qw/F R U/); 
				print {$outFh} join("\t", @{$ta[$i]}[0,1,2,3,4], @{$ta[$j]}[0,1,2,3,4], $blk_num, $str_i, $str_j)."\n"; 
				
			} else {
				print {$outFh} join("\t", @{$ta[$i]}[0,1,2,3,4], @{$ta[$j]}[0,1,2,3,4], $blk_num)."\n"; 
			}
		}
	}


	return; 
}# cnvt_coords_gff2tab() 

sub fix_strand_gff {
	my %data = %{ &_load_coords_gff( $opts{'in_blkGff'} ) }; 
	my %scfID_order = %{ &_get_scfID_order( \%data ) }; 
	&_fix_blk_strand( \%data, \%scfID_order ); 
	&_out_coords_gff( \%data, $opts{'out_blk'} ); 
}# fix_strand_gff() 

sub sort_coords_gff {
	my %data = %{ &_load_coords_gff( $opts{'in_blkGff'} ) }; 
	my %scfID_order = %{ &_get_scfID_order( \%data ) }; 
	my @sorted_blk_num = @{ &_sort_blkNum_by_IdPos( \%data, \%scfID_order ) }; 
	&_out_coords_gff( \%data, $opts{'out_blk'}, \@sorted_blk_num ); 
	return; 
}# sort_coords_gff() 
sub cnvt_coords_txt2gff {
	my %data = %{ &_load_coords_txt( $opts{'in_blkTxt'} ) }; 
	&_out_coords_gff( \%data , $opts{'out_blk'} ); 
	return; 
}# cnvt_coords_txt2gff () 
sub cnvt_coords_gff2txt {
	my %data = %{ &_load_coords_gff( $opts{'in_blkGff'} ) }; 
	&_out_coords_txt( \%data , $opts{'out_blk'} ); 
}# cnvt_coords_gff2txt ()

################################################################################
# Inner sub-routines. 
################################################################################

=head1 _out_coords_txt ( \%coords_hash, [*STDOUT | $out_fh], [@blk_num_to_output] ) 
=cut
sub _out_coords_txt {
	my ($href, $outFh, $blk_num_aref) = @_; 
	$outFh //= \*STDOUT; 
	( defined $blk_num_aref and @{$blk_num_aref} > 0 ) or $blk_num_aref = [ sort { $href->{'blk_num2order'}{$a} <=> $href->{'blk_num2order'}{$b} } keys %{$href->{'blk_num2order'}} ]; 
	print {$outFh} join("\t", qw/Seq_id  Size    Description/)."\n"; 
	for my $scfID (sort { $href->{'scfID2seq_id'}{$a} <=> $href->{'scfID2seq_id'}{$b} } keys %{$href->{'scfID2seq_id'}}) {
		print {$outFh} join("\t", $href->{'scfID2seq_id'}{$scfID}, $href->{'scfID2len'}{$scfID}, $scfID)."\n"; 
	}
	print {$outFh} '-' x 80 . "\n"; 
	for my $blk_num (@$blk_num_aref) {
		print {$outFh} "Block #$blk_num\n"; 
		print {$outFh} join("\t", qw/Seq_id  Strand  Start   End     Length/)."\n"; 
		for my $a1 ( @{$href->{'blk_num2inf'}{$blk_num}} ) {
			my @ta = @$a1; 
			$ta[1] eq '-' and do { @ta[2,3]=@ta[3,2]; $ta[2] >= $ta[3] or &stopErr("[Err] Bad format of input info [@ta[2,3,1]]\n"); }; 
			my $seq_id = $href->{'scfID2seq_id'}{$ta[0]}; 
			print {$outFh} join("\t", $seq_id, @ta[1,2,3,4])."\n"; 
		}
		print {$outFh} '-' x 80 . "\n"; 
	}
	return; 
}# _out_coords_txt () 

=head1 _out_coords_gff ( \%coords_hash, [*STDOUT | $out_fh], [@blk_num_to_output] )

Format of blocks_coords.gff : 
 ##gff-version 2
 ##source-version Sibelia 3.0.6
 ##Type DNA
 Cmo_Chr03       Sibelia synteny_block_copy      2109420 2533163 .       -       .       1
 Cma_Chr03       Sibelia synteny_block_copy      1341335 1727105 .       +       .       1
 Cma_Chr03       Sibelia synteny_block_copy      2170358 2338506 .       +       .       2
 Cmo_Chr03       Sibelia synteny_block_copy      3214481 3364918 .       +       .       2
 Cmo_Chr03       Sibelia synteny_block_copy      2573819 2579977 .       +       .       3
 Cma_Chr03       Sibelia synteny_block_copy      1296357 1302498 .       -       .       3
 ... 

=cut
sub _out_coords_gff {
	my ($href, $outFh, $blk_num_aref) = @_; 
	$outFh //= \*STDOUT; 
	( defined $blk_num_aref and @{$blk_num_aref} > 0 ) or $blk_num_aref = [ sort { $href->{'blk_num2order'}{$a} <=> $href->{'blk_num2order'}{$b} } keys %{$href->{'blk_num2order'}} ]; 
	for my $blk_num (@$blk_num_aref) {
		my @ta = @{$href->{'blk_num2inf'}{$blk_num}}; 
		for my $a1 ( @{$href->{'blk_num2inf'}{$blk_num}} ) {
			my @ta = @$a1; 
			my $ele_ID = "$ta[0]_${blk_num}"; 
			if (defined $glob{'has_ID'}{$ele_ID}) {
				my $copyN = 1; 
				my $new_ele_ID = sprintf("%s_dup%03d", "${ele_ID}", "$copyN"); 
				while ( defined $glob{'has_ID'}{$new_ele_ID} ) {
					$copyN ++; 
					$copyN > 1e3 and &stopErr("[Err] Failed to find an ID for [$ele_ID]\n"); 
					$new_ele_ID = sprintf("%s_dup%03d", "${ele_ID}", "$copyN"); 
				}
				$ele_ID = $new_ele_ID; 
			}
			$glob{'has_ID'}{$ele_ID} = 1; 
			my $len = $ta[3]-$ta[2]+1; 
			print {$outFh} join("\t", $ta[0], 'Sibelia', 'synteny_block_copy', @ta[2,3], '.', $ta[1], '.', "ID=${ele_ID};block_num=$blk_num;ele_len=$len")."\n"; 
		}
	}
	return; 
}# _out_coords_gff() 

=head1 _load_coords_txt ( 'blocks_coords.txt' ) 

Format of blocks_coords.txt : 
 Seq_id  Size    Description
 1       35366548        WM97v2_Chr01
 2       35446621        WMPI_Chr01
 --------------------------------------------------------------------------------
 Block #1
 Seq_id  Strand  Start   End     Length
 1       +       12120463        12226801        106339
 2       -       13521828        13434826        87003
 --------------------------------------------------------------------------------
 Block #2
 Seq_id  Strand  Start   End     Length
 1       -       11167513        11044299        123215
 2       +       10847143        10969522        122380
 --------------------------------------------------------------------------------
 .....

Return    : (\%info)
  $info{'blk_num2inf'}   = {$blk_num} => [[ $scfID, $strand(+/-), $start, $end, $length ], [], ...]
    Here when $strand eq '-', $start should be bigger than $end, and they are 1-based indexed. 
  $info{'scfID2len'}     = {$scfID} => $scaffold_length;
  $info{'scfID2seq_id'}  = seq_id;
  $info{'blk_num2order'} = sequence_order_of_blk_num;

=cut
sub _load_coords_txt {
	my $fn = shift; 
	my $fh = &_openFnOrSTDIO( $fn, '<' ); 
	my %back; 
	my %tag; 
	$tag{'is_blk'} = 0; 
	while (<$fh>) {
		chomp; 
		m/^\-*$/ and next; 
		my @ta = split(/\t/, $_); 
		$ta[0] eq 'Seq_id' and next; 
		if ( $tag{'is_blk'} and @ta >= 5 ) {
			my $scfID = $tag{'num2inf'}{$ta[0]}[1]; 
			defined $scfID or &stopErr("[Err] Undefined seq_id [$ta[0]]\n"); 
			$ta[1] eq '-' and @ta[2,3] = @ta[3,2]; 
			$ta[2] <= $ta[3] or &stopErr("[Err] Bad number order: [@ta[2,3,1]]\n"); 
			push( @{$back{'blk_num2inf'}{$tag{'blk_num'}}}, [ $scfID, @ta[1 .. 4] ] ); 
		} elsif ( @ta == 3 ) {
			$ta[0] =~ m!^\d+$! or &stopErr("[Err] Bad seq_id [$ta[0]]\n"); 
			$tag{'num2inf'}{$ta[0]} = [$ta[1], $ta[2]]; # { id_num } => [ seq_size, seq_ID ]
			$back{'scfID2len'}{$ta[2]} = $ta[1]; 
			$back{'scfID2seq_id'}{$ta[2]} = $ta[1]; 
			defined $tag{'has_scfID'}{$ta[2]} and &stopErr("[Err] repeated scaffold ID [$ta[2]]\n"); 
			$tag{'has_scfID'}{$ta[2]} = 1; 
		} elsif ( m/^Block\s+\#(\d+)\s*$/ ) {
			$tag{'blk_num'} = $1; 
			$tag{'is_blk'} = 1; 
			$tag{'blk_num_order'} ++; 
			$back{'blk_num2order'}{$tag{'blk_num'}} = $tag{'blk_num_order'}; 
		} else {
			&stopErr("[Err] Bad line: $_\n"); 
		}
	}
	close($fh); 


	return(\%back); 
}# _load_coords_txt() 

=head1 _load_coords_gff ( 'blocks_coords.gff' ) 

Return same hash reference with _load_coords_txt(), but read in .gff file. 

=cut
sub _load_coords_gff {
	my $fn = shift; 
	my $fh = &_openFnOrSTDIO( $fn, '<' ); 
	my %back; 
	my %cnt; 
	while (&wantLineC($fh)) {
		my @ta = split("\t", $_); 
		$back{'scfID2len'}{$ta[0]} //= 'NA'; 
		unless ( defined $back{'scfID2seq_id'}{$ta[0]} ) {
			$cnt{'seq_id'} ++; 
			$back{'scfID2seq_id'}{$ta[0]} = $cnt{'seq_id'}; 
		}
		my $blk_num; 
		if ($ta[8] =~ m!^(\d+)$!) {
			$blk_num = $1; 
		} elsif ( $ta[8] =~ m!(?:^\s*|;\s*)ID=(\S+)_(\d+)(?:\s*$|\s*;)!i ) {
			$ta[0] eq $1 or &stopErr("[Err] Bad line: $_\n"); 
			$blk_num = $2; 
		} elsif ( $ta[8] =~ m!(?:^\s*|;\s*)ID=(\S+)_(\d+)_dup(\d+)(?:\s*$|\s*;)!i ) {
			$ta[0] eq $1 or &stopErr("[Err] Bad line: $_\n"); 
			$blk_num = $2; 
		} else {
			&stopErr("[Err] Failed to parse line: $_\n"); 
		}
		$blk_num =~ s!^0+!!g; 
		push(@{$back{'blk_num2inf'}{$blk_num}}, [ $ta[0], $ta[6], $ta[3], $ta[4], $ta[4]-$ta[3]+1 ]); 
		unless ( defined $back{'blk_num2order'}{$blk_num} ) {
			$cnt{'blk_num'} ++; 
			$back{'blk_num2order'}{$blk_num} = $cnt{'blk_num'}; 
		}
	}
	close($fh); 
	&_add_scfID2len( \%back ); 
	return(\%back); 
}# _load_coords_gff () 

sub _openFnOrSTDIO {
	my ($fn, $type) = @_; 
	$type //= '<'; 
	if (defined $fn and $fn ne '') {
		return(&openFH($fn, $type)); 
	} elsif ( $type eq '>' ) {
		return( \*STDOUT ); 
	} elsif ( @ARGV > 0 ) {
		$fn = shift(@ARGV); 
		return(&openFH($fn, $type)); 
	} elsif ( !(-t) ) {
		if ( $type eq '<' ) {
			return( \*STDIN ); 
		} else {
			&stopErr("[Err] Bad type [$type]\n"); 
		}
	}
	&stopErr("[Err] Failed to parse ($fn, $type)\n"); 
	return; 
}# _openFnOrSTDIO() 

sub _add_scfID_order {
	my ($href, $scfID) = @_; 
	defined $href->{$scfID} and return; 
	my $n = scalar(keys %$href); 
	$href->{$scfID} = $n+1; 
	return; 
}# _add_scfID_order ()

=head1 _get_scfID_order( \%info )

 Function   : Use -scfID_order_list , -scfID_order@ and %info to get order of each scfID;
 Return     : (\%scfID_order)
 $scfID_order{$scfID} = $order_number ; # $order_number is from 1 to the end.

=cut
sub _get_scfID_order {
	my ($href) = @_; 
	defined $href or &stopErr("[Err] No input for \&_get_scfID_order()\n"); 
	my %back; 

	if ( defined $opts{'scfID_order_list'} ) {
		my $fh = &openFH($opts{'scfID_order_list'}, '<'); 
		while (&wantLineC($fh)) {
			my @ta = &splitL("\t", $_); 
			my $id = $ta[0]; 
			&_add_scfID_order( \%back, $id ); 
		}
		close($fh); 
	}
	if ( defined $opts{'scfID_order'} ) { 
		for my $id (@{$opts{'scfID_order'}}) {
			&_add_scfID_order( \%back, $id); 
		}
	}
	for my $id ( grep { ! defined $back{$_} } sort { $href->{'scfID2seq_id'}{$a} <=> $href->{'scfID2seq_id'}{$b} } keys %{ $href->{'scfID2seq_id'} } ) {
		&_add_scfID_order( \%back, $id ); 
	}

	return(\%back); 
}# _get_scfID_order() 

=head1 _blkNum_represent_val( \%info , \%scfID_order )

Return      : (\%blkNum_to_val)
 $blkNum_to_val{$blk_num} = [ $first_scfID, $first_scf_posi_start, $first_scf_posi_end ]

=cut
sub _blkNum_represent_val {
	my ($href_blk, $href_order) = @_; 
	$href_order //= &_get_scfID_order( $href_blk ); 
	my %val; 
	for my $blk_num (sort { $href_blk->{'blk_num2order'}{$a} <=> $href_blk->{'blk_num2order'}{$b} } keys %{$href_blk->{'blk_num2order'}}) {
		my $ref_scfID = ( sort { $href_order->{$a} <=> $href_order->{$b} } map { $_->[0] } @{$href_blk->{'blk_num2inf'}{$blk_num}} )[0]; 
		$val{$blk_num} = [ $ref_scfID ]; 
		for my $a1 ( @{$href_blk->{'blk_num2inf'}{$blk_num}} ) {
			$a1->[0] eq $ref_scfID or next; 
			$val{$blk_num}[1] //= $a1->[2]; 
			$val{$blk_num}[2] //= $a1->[3]; 
			if ( $val{$blk_num}[1] > $a1->[2] ) {
				$val{$blk_num}[1] = $a1->[2]; 
				$val{$blk_num}[2] = $a1->[3]; 
			} elsif ( $val{$blk_num}[1] == $a1->[2] and $val{$blk_num}[2] > $a1->[3] ) {
				$val{$blk_num}[1] = $a1->[2]; 
				$val{$blk_num}[2] = $a1->[3]; 
			}
		}
	}

	return(\%val); 
}# _blkNum_represent_val () 

=head1 _fix_pairScf_order( \%info, \%scfID_order )

Function   : If the second seq_id is in front of its paired segment, those two will be replaced. 

=cut
sub _fix_pairScf_order {
	my ($href_blk, $href_order) = @_; 
	$href_order //= &_get_scfID_order( $href_blk ); 

	for my $blk_num ( sort { $href_blk->{'blk_num2order'}{$a} <=> $href_blk->{'blk_num2order'}{$b} } keys %{$href_blk->{'blk_num2order'}} ) {
		my @ta = @{$href_blk->{'blk_num2inf'}{$blk_num}}; 
		for (my $i=0; $i<@ta; $i+=2) {
			my $j=$i+1; 
			defined $ta[$j] or last; 
			@ta[$i,$j] = sort { $href_order->{$a->[0]} <=> $href_order->{$b->[0]} || $a->[2] <=> $b->[2] || $a->[3] <=> $a->[3] } @ta[$i,$j]; 
		}
		@{$href_blk->{'blk_num2inf'}{$blk_num}} = @ta; 
	}
	return; 
}# _fix_pairScf_order() 

=head1 _fix_blk_strand( \%info, \%scfID_order )

Function    : If the representative segment of a block is '-' strand, all segments in this block will be reversed. 

Return      : () 

=cut
sub _fix_blk_strand {
	my ($href_blk, $href_order) = @_; 
	$href_order //= &_get_scfID_order( $href_blk ); 

	for my $blk_num ( sort { $href_blk->{'blk_num2order'}{$a} <=> $href_blk->{'blk_num2order'}{$b} } keys %{$href_blk->{'blk_num2order'}} ) {
		my $ref_segment_aref = ( sort { $href_order->{$a->[0]} <=> $href_order->{$b->[0]} || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3] || $a->[1] cmp $b->[1] } @{$href_blk->{'blk_num2inf'}{$blk_num}} )[0]; 
		if ( $ref_segment_aref->[1] eq '-' ) {
			for my $a1 ( @{$href_blk->{'blk_num2inf'}{$blk_num}} ) {
				$a1->[1] =~ tr/+-/-+/; 
			}
		}
	}
	return; 
}# _fix_blk_strand() 

=head1 _sort_blkNum_by_IdPos( \%info, \%scfID_order )

Return      : ( \@sorted_blk_num )
 block numbers (blk_num) are sorted by [$scfID, $startPos, $endPos]. 

=cut
sub _sort_blkNum_by_IdPos {
	my ($href_blk, $href_order) = @_; 
	my @back; 

	my %val = %{ &_blkNum_represent_val( $href_blk, $href_order ) }; 

	@back = sort { $href_order->{$val{$a}[0]} <=> $href_order->{$val{$b}[0]} || $val{$a}[1] <=> $val{$b}[1] || $val{$a}[2] <=> $val{$b}[2] } keys %val; 

	return(\@back); 
}# _sort_blkNum_by_IdPos ()

=head1 _add_scfID2len( \%blk_hash )
=cut
sub _add_scfID2len {
	my ($href_blk) = @_; 
	defined $opts{'in_scfLen'} or return; 
	my %scfID2len = %{ &_load_keyLen( $opts{'in_scfLen'} ) }; 
	for my $scfID (keys %{$href_blk->{'scfID2len'}}) {
		$href_blk->{'scfID2len'}{$scfID} eq 'NA' and defined $scfID2len{$scfID} and $href_blk->{'scfID2len'}{$scfID} = $scfID2len{$scfID}; 
	}
	return; 
}# scfID2len() 

=head1 _load_keyLen( 'in.fa.key_len' )
=cut
sub _load_keyLen {
	my $fn = shift; 
	my $fh = &openFH($fn, '<'); 
	my %back; 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		$back{$ta[0]} //= $ta[1]; 
	}
	close($fh); 
	return(\%back); 
}# _load_keyLen() 

