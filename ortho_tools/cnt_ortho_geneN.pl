#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"in_file:s", # all_orthomcl.out
	"in_fmt:s", # Default orthomcl
	"c1_required_taxaLis:s@", "c1_required_taxaNum:i@", 
	"c2_min_taxaN:i", # At least 
	"c3_chk_taxaLis:s", "c3_min_taxaN:i", 
	"c3_bad_taxaLis:s", "c3_max_taxaN:i", 
	"c4_gg_list:s", # all_12_taxa.gg 
	"c5_inner_taxaLis:s", "c5_outer_taxaLis:s", 
	"c6_inner_taxaLis:s", "c6_outer_taxaLis:s", 
	"c9_ref_taxaLis:s@", # 
	"odir:s", 
); 

my $help_txt = <<HH; 
################################################################################
# perl $0 -in_file all_orthomcl.out -odir out_new_dir
#
# -in_fmt                 ['orthomcl'] Could be : orthomcl , orthofinder ; 
# 
# category_01 : All lists must be satisfied. 
#  -c1_required_taxaLis   [taxaList\@] One taxID per line. 
#  -c1_required_taxaNum   [taxaNum\@] Number of taxa needed in related list
#
# category_02 : 
#  -c2_min_taxaN          [taxaNum] Number of taxa needed in a OG; 
#
# category_03 : 
#  -c3_chk_taxaLis        [taxaList] Subset of taxa to check 
#  -c3_min_taxaN          [taxaNum] Number of taxa needed in subset. 
#  -c3_bad_taxaLis        [taxaList] Subset of taxa to keep out
#  -c3_max_taxaN          [taxaNum] Number of taxa allowed in subset. 
#
# category_04 : 
#  -c4_gg_list            [all_12_taxa.gg] Format from OrthoMCL; 
#
# category_05 : 
#  -c5_inner_taxaLis      [taxaList] Wanted taxa list; 
#  -c5_outer_taxaLis      [taxaList] Excluded taxa list; 
#
# category_06 : 
#  -c6_inner_taxaLis      [taxaList] Wanted taxa list; 
#  -c6_outer_taxaLis      [taxaList] Excluded taxa list; 
#
# category_07 : 
#
# category_08 : 
#
# category_09 : At least one list must be satisfied. 
#  -c9_ref_taxaLis        [taxaList\@] One taxID per line. 
# 
################################################################################
# category_01 : At least \$c1_required_taxaNum in \%c1_required_taxaLis; 
# category_02 : Not in previous categories, at least \$c2_min_taxaN 
# category_03 : Not in previous categories, at least \$c3_min_taxaN taxa within \%c3_chk_taxaLis , and at most \$c3_max_taxaN within \%c3_bad_taxaLis; 
# category_04 : Not in previous categories, without any homolog in other species; 
# category_05 : Not in previous categories, have all in \%c5_inner_taxaLis and non in \%c5_outer_taxaLis; 
# category_06 : Not in previous categories, have all in \%c6_inner_taxaLis and non in \%c6_outer_taxaLis; 
# category_07 : Not in previous categories; 
# category_08 : At least one paralog within species; 
# category_09 : Having all homolog with at least one taxaLis in \@c8_ref_taxaLis; 

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
defined $opts{'in_file'} or &LogInforSunhh::usage($help_txt); 
$opts{'odir'} //= '.'; 

my %glob; 
&prepare_glob(); 

my ( @ortho_grps , %taxGen2grpID); 
{
	my ($r1, $r2) = &load_grps( $opts{'in_file'} , $glob{'in_fmt'} ); 
	@ortho_grps = @$r1; 
	%taxGen2grpID = %$r2; 
}

&cnt_c1();  
&cnt_c2(); 
&cnt_c3(); 
&cnt_c4(); 
&cnt_c5(); 
&cnt_c6(); 
&cnt_c7(); 
&cnt_c8(); 
&cnt_c9(); 
&out_cN_geneN([qw/1 2 3 4 5 6 7 8 9/]); 

=head1 cnt_c1() 

category_01 : At least $c1_required_taxaNum in %c1_required_taxaLis; 
 Return : ( )
 Modify : %glob keys qw/c1_required_taxaLis c1_required_taxaNum c1_required_taxaAll c1_has_taxGene c1_has_grpID c1_has_geneN/
   $glob{'c1_required_taxaLis'} = [ { taxa_ID => $lineN }, {}, ... ]; 
   $glob{'c1_required_taxaNum'} = [ $opts{'c1_required_taxaNum'}[0], $opts{'c1_required_taxaNum'}[1], ... ]; 
   $glob{'c1_required_taxaAll'} = [ $total_required_taxaN[0], $total_required_taxaN[1], ... ];
   update %glob keys qw/c${cN}_has_taxGene c${cN}_has_grpID c${cN}_has_geneN/

=cut
sub cnt_c1 {
	$glob{'c1_required_taxaLis'} //= [{}]; 
	$glob{'c1_required_taxaNum'} //= [0]; 
	if ( defined $opts{'c1_required_taxaLis'} ) {
		for (my $i=0; $i<@{$opts{'c1_required_taxaLis'}}; $i++) {
			( $glob{'c1_required_taxaLis'}[$i] ) = &load_taxaLis( $opts{'c1_required_taxaLis'}[$i] ); 
			$glob{'c1_required_taxaAll'}[$i] = scalar(keys %{$glob{'c1_required_taxaLis'}[$i]}); 
			$opts{'c1_required_taxaNum'}[$i] //= $glob{'c1_required_taxaAll'}[$i]; 
			$glob{'c1_required_taxaNum'}[$i] = $opts{'c1_required_taxaNum'}[$i]; 
		}
	}
	my %cnt; 
	for my $r1 (@ortho_grps) {
		my @taxa_ids = keys %{$r1->{'tax2gene'}}; 
		my $is_good = 1; 
		for ( my $i=0; $i< @{$glob{'c1_required_taxaLis'}}; $i++ ) {
			$cnt{'curr_requiredN'} = 0; 
			for my $t1 (@taxa_ids) {
				defined $glob{'c1_required_taxaLis'}[$i]{$t1} and $cnt{'curr_requiredN'}++; 
			}
			$cnt{'curr_requiredN'} >= $glob{'c1_required_taxaNum'}[$i] or do { $is_good = 0; last; }; 
		}
		if ( $is_good == 1 ) {
			&fill_has( '1', $r1 ); 
		}
	}

	&out_data('1'); 
	return; 
}# cnt_c1 ()

=head1 cnt_c2() 

category_02 : Not in previous categories, at least $c2_min_taxaN 

=cut
sub cnt_c2 {
	defined $opts{'c2_min_taxaN'} or do { $opts{'c2_min_taxaN'} = 1; &tsmsg("[Wrn] set default -c2_min_taxaN to $opts{'c2_min_taxaN'}\n");  }; 
	for my $r1 (@ortho_grps) {
		&chk_prev_has( 'grpID', ['1'], $r1->{'grpID'} ) and next; 
		my $is_good = 1; 
		my $tn = scalar( keys %{$r1->{'tax2gene'}} ); 
		$tn >= $opts{'c2_min_taxaN'} or $is_good = 0; 
		$is_good == 1 or next; 
		&fill_has( '2', $r1 ); 
	}

	&out_data('2'); 
	return; 
}# cnt_c2 () 



=head1 cnt_c3()

category_03 : Not in previous categories, at least $c3_min_taxaN taxa within %c3_chk_taxaLis ; 

=cut
sub cnt_c3 {
	$glob{'c3_chk_taxaLis'} //= {}; 
	$glob{'c3_min_taxaN'}   //= 0; 
	$glob{'c3_bad_taxaLis'} //= {}; 
	$glob{'c3_max_taxaN'}   //= 0; 
	if ( defined $opts{'c3_chk_taxaLis'} ) {
		( $glob{'c3_chk_taxaLis'} ) = &load_taxaLis( $opts{'c3_chk_taxaLis'} ); 
		$glob{'c3_chk_taxaAll'} = scalar( keys %{$glob{'c3_chk_taxaLis'}} ); 
		$opts{'c3_min_taxaN'} //= $glob{'c3_chk_taxaAll'}; 
		$glob{'c3_min_taxaN'} = $opts{'c3_min_taxaN'}; 
	}
	if ( defined $opts{'c3_bad_taxaLis'} ) {
		( $glob{'c3_bad_taxaLis'} ) = &load_taxaLis( $opts{'c3_bad_taxaLis'} ); 
		# $glob{'c3_bad_taxaAll'} = scalar( keys %{$glob{'c3_bad_taxaLis'}} ); 
		$glob{'c3_max_taxaN'} = $opts{'c3_max_taxaN'}; 
	}

	for my $r1 (@ortho_grps) {
		&chk_prev_has( 'grpID', ['1', '2'], $r1->{'grpID'} ) and next; 
		my @taxa_ids = keys %{ $r1->{'tax2gene'} }; 
		my $is_good = 1; 
		my %cnt = ( 'curr_requiredN' => 0 , 'curr_badN' => 0 ); 
		for my $t1 (@taxa_ids) {
			defined $glob{'c3_chk_taxaLis'}{$t1} and $cnt{'curr_requiredN'}++; 
			defined $glob{'c3_bad_taxaLis'}{$t1} and $cnt{'curr_badN'}++; 
		}
		$cnt{'curr_badN'} <= $glob{'c3_max_taxaN'} or $is_good = 0; 
		$cnt{'curr_requiredN'} >= $glob{'c3_min_taxaN'} or $is_good = 0; 
		$is_good == 1 or next; 
		&fill_has( '3', $r1 ); 
	}

	&out_data('3'); 
	return; 
}#cnt_c3 () 

=head1 cnt_c6()

category_06 : Not in previous categories, have all in %c6_inner_taxaLis and non in %c6_outer_taxaLis; 

=cut
sub cnt_c6 {
	$glob{'c6_inner_taxaLis'} = {}; 
	$glob{'c6_inner_taxaAll'} = 0; 
	$glob{'c6_outer_taxaLis'} = {}; 
	if ( defined $opts{'c6_inner_taxaLis'} ) {
		$glob{'c6_inner_taxaLis'} = &load_taxaLis( $opts{'c6_inner_taxaLis'} ); 
		$glob{'c6_inner_taxaAll'} = scalar( keys %{$glob{'c6_inner_taxaLis'}} ); 
	}
	defined $opts{'c6_outer_taxaLis'} and $glob{'c6_outer_taxaLis'} = &load_taxaLis( $opts{'c6_outer_taxaLis'} ); 

	for my $r1 (@ortho_grps) {
		&chk_prev_has( 'grpID', ['1', '2', '3', '4','5'], $r1->{'grpID'} ) and next; 
		my @taxa_ids = keys %{ $r1->{'tax2gene'} }; 
		my $is_good = 1; 
		my %cnt = ( 'curr_requiredN' => 0 ); 
		for my $t1 (@taxa_ids) {
			defined $glob{'c6_outer_taxaLis'}{$t1} and do { $is_good = 0; last; }; 
			defined $glob{'c6_inner_taxaLis'}{$t1} and $cnt{'curr_requiredN'}++; 
		}
		$cnt{'curr_requiredN'} == $glob{'c6_inner_taxaAll'} or $is_good = 0; 
		$is_good == 1 or next; 
		&fill_has( '6', $r1 );
	}
	&out_data('6'); 
	return; 
}# cnt_c6 () 

=head1 cnt_c5()

category_05 : Not in previous categories, have all in %c5_inner_taxaLis and non in %c5_outer_taxaLis; 

=cut
sub cnt_c5 {
	$glob{'c5_inner_taxaLis'} = {}; 
	$glob{'c5_inner_taxaAll'} = 0; 
	$glob{'c5_outer_taxaLis'} = {}; 
	if ( defined $opts{'c5_inner_taxaLis'} ) {
		$glob{'c5_inner_taxaLis'} = &load_taxaLis( $opts{'c5_inner_taxaLis'} ); 
		$glob{'c5_inner_taxaAll'} = scalar( keys %{$glob{'c5_inner_taxaLis'}} ); 
	}
	defined $opts{'c5_outer_taxaLis'} and $glob{'c5_outer_taxaLis'} = &load_taxaLis( $opts{'c5_outer_taxaLis'} ); 

	for my $r1 (@ortho_grps) {
		&chk_prev_has( 'grpID', ['1', '2', '3', '4'], $r1->{'grpID'} ) and next; 
		my @taxa_ids = keys %{ $r1->{'tax2gene'} }; 
		my $is_good = 1; 
		my %cnt = ( 'curr_requiredN' => 0 ); 
		for my $t1 (@taxa_ids) {
			defined $glob{'c5_outer_taxaLis'}{$t1} and do { $is_good = 0; last; }; 
			defined $glob{'c5_inner_taxaLis'}{$t1} and $cnt{'curr_requiredN'}++; 
		}
		$cnt{'curr_requiredN'} == $glob{'c5_inner_taxaAll'} or $is_good = 0; 
		$is_good == 1 or next; 
		&fill_has( '5', $r1 );
	}
	&out_data('5'); 
	return; 
}# cnt_c5 () 

=head1 cnt_c4()

category_04 : Not in previous categories, without any homolog in other species; 
              Including genes not captured in ortholog groups. 

=cut
sub cnt_c4 {
	$glob{'c4_gg_list'} = {}; 
	if ( defined $opts{'c4_gg_list'} ) {
		$glob{'c4_gg_list'} = &load_ggList( $opts{'c4_gg_list'} ); 
	}
	my %processed; 
	for my $r1 (@ortho_grps) {
		my @taxa_ids = keys %{ $r1->{'tax2gene'} }; 
		for my $t1 (@taxa_ids) {
			for my $t2 (keys %{$r1->{'tax2gene'}{$t1}}) {
				$processed{$t1}{$t2} = 1; 
			}
		}
		&chk_prev_has( 'grpID', ['1','2','3'], $r1->{'grpID'} ) and next; 
		$#taxa_ids > 0 and next; 
		&fill_has( '4', $r1 ); 
	}

	for my $t1 ( keys %{$glob{'c4_gg_list'}{'tax2gene_a'}} ) {
		for my $t2 ( @{ $glob{'c4_gg_list'}{'tax2gene_a'}{$t1} } ) {
			defined $processed{$t1}{$t2} and next; 
			$glob{'c4_has_geneN'}{$t1} ++; 
			$glob{'c4_has_taxGene'}{$t1}{$t2} = -1; 
		}
	}

	&out_data('4'); 
	return; 
}# cnt_c4() 

=head1 cnt_c7()

category_07 : Not in previous categories; 

=cut
sub cnt_c7 {
	for my $r1 ( @ortho_grps ) {
		&chk_prev_has( 'grpID', [qw/1 2 3 4 5 6/], $r1->{'grpID'} ) and next; 
		for my $t1 (keys %{$r1->{'tax2gene'}}) {
			for my $t2 ( keys %{$r1->{'tax2gene'}{$t1}} ) {
				&chk_prev_has( 'taxGene', ['5'], [ $t1, $t2 ] ) and &stopErr("[Err] Bad category for gene {$t1}{$t2}\n"); 
				&chk_prev_has( 'taxGene', ['6'], [ $t1, $t2 ] ) and &stopErr("[Err] Bad category for gene {$t1}{$t2}\n"); 
			}
		}
		&fill_has( '7', $r1 ); 
	}
	&out_data('7'); 
	return; 
}# cnt_c7 () 

=head1 cnt_c8 () 

category_08 : At least one paralog within species; 

=cut
sub cnt_c8 {
	for my $r1 (@ortho_grps) {
		my @taxa_ids = keys %{$r1->{'tax2gene'}}; 
		for my $t1 (@taxa_ids) {
			my @gene_ids = keys %{$r1->{'tax2gene'}{$t1}} ; 
			scalar(@gene_ids) > 1 or next; 
			$glob{'c8_has_grpID'}{$r1->{'grpID'}} = 1; 
			for my $t2 (@gene_ids) {
				$glob{'c8_has_geneN'}{$t1} ++; 
				$glob{'c8_has_taxGene'}{$t1}{$t2} = $r1->{'tax2gene'}{$t1}{$t2}; 
			}
		}
	}
	&out_data('8'); 
	return; 
}# cnt_c8 () 

=head1 cnt_c9 () 

category_09 : Having homolog with at least one taxaLis in @c9_ref_taxaLis; 

=cut
sub cnt_c9 {
	$glob{'c9_ref_taxaLis'} = []; 
	$glob{'c9_ref_taxaAll'} = [0];  
	if ( defined $opts{'c9_ref_taxaLis'} ) {
		for (my $i=0; $i<@{$opts{'c9_ref_taxaLis'}}; $i++) {
			$glob{'c9_ref_taxaLis'}[$i] = &load_taxaLis( $opts{'c9_ref_taxaLis'}[$i] ); 
			$glob{'c9_ref_taxaAll'}[$i] = scalar( keys %{$glob{'c9_ref_taxaLis'}[$i]} ); 
		}
	}
	for my $r1 (@ortho_grps) {
		my @taxa_ids = keys %{$r1->{'tax2gene'}}; 
		my $is_good = 0; 
		for ( my $i=0; $i<@{$glob{'c9_ref_taxaLis'}}; $i++ ) {
			my %cnt = ( 'curr_requiredN' => 0 ); 
			for my $t1 ( @taxa_ids ) {
				defined $glob{'c9_ref_taxaLis'}[$i]{$t1} and $cnt{'curr_requiredN'}++; 
			}
			$cnt{'curr_requiredN'} == $glob{'c9_ref_taxaAll'}[$i] and do { $is_good = 1; last; }; 
		}
		$is_good == 1 or next; 
		&fill_has( '9', $r1 ); 
	}
	&out_data('9'); 
	return; 
}# cnt_c9 () 

=head1 fill_has ( $cN, $ortho_grps[$selected] )

Function : Fill $glob{"c${cN}_has_taxGene"} and $glob{"c${cN}_has_grpID"} and $glob{"c${cN}_has_geneN"}

Return   : 
   $glob{"c${cN}_has_taxGene"}      = { taxa_ID => { gene_ID => $grp_order , gene_ID=>$grp_order, ... }, ... }; 
   $glob{"c${cN}_has_grpID"}        = { $grpID  => 1 }
   $glob{"c${cN}_has_geneN"}        = { taxa_ID => $number_of_c1_genes }

=cut
sub fill_has {
	my ($cN, $og) = @_; 
	$glob{"c${cN}_has_grpID"}{ $og->{'grpID'} } = 1; 
	my $k1 = "c${cN}_has_taxGene"; 
	my $k2 = "c${cN}_has_geneN"; 
	for my $t1 ( keys %{$og->{'tax2gene'}} ) {
		for my $t2 ( keys %{$og->{'tax2gene'}{$t1}} ) {
			$glob{$k1}{$t1}{$t2} = $og->{'tax2gene'}{$t1}{$t2}; 
			$glob{$k2}{$t1} ++; 
		}
	}

	return; 
}# fill_has ()

=head1 load_grps( $opts{'in_file'} )

Return : ( [ \%grp_1, \%grp_2, ... ], %gene2grpID )
  keys %grp = qw/grpID geneN taxaN tax2gene/
    $grp{'grpID'} = ORTHOMCL\d+ 
    $grp{'geneN'} = \d+ 
    $grp{'taxaN'} = \d+
    $grp{'tax2gene'} = { {taxa_1}{gene_1} => $order, {taxa_1}{gene_2} => $order, {taxa_2}{gene_3} => $order, ... }
  keys %gene2grpID = ( gene_1, gene_2, gene_3, ... ) 
    $gene2grpID{$taxa_1}{$gene_1} = $grpID

=cut
sub load_grps {
	my ($fn, $fmt) = @_; 
	$fmt //= 'orthocml'; 
	if ($fmt =~ m!^\s*orthomcl\s*$!i) {
		return( &load_grps_fmt_orthocml($fn) ); 
	} elsif ($fmt =~ m!^\s*orthofinder\s*$!i) {
		return( &load_grps_fmt_orthofinder($fn) ); 
	} else {
		&stopErr("[Err] Unknown format of input Orthologous Group file.\n"); 
	}
	return; 
}# load_grps() 

sub load_grps_fmt_orthofinder {
	my $fn = shift; 
	my $fh = &openFH($fn, '<'); 
	my (@back, %gene2grpID); 
	# For OrthologousGroups.csv 
	my @hh; 
	while (&wantLineC($fh)) {
		chomp; 
		my @ta = &splitL("\t", $_); 
		if ($. == 1) {
			$ta[0] eq '' or &stopErr("[Err] Bad 1st line: $_\n"); 
			@hh = @ta; 
			next; 
		}
		my ($grpID, $geneN, $taxaN) = ($ta[0], 0, 0); 
		push(@back, {}); 
		$back[-1]{'grpID'} = $grpID; 
		for (my $i=1; $i<@ta; $i++) {
			$ta[$i] eq '' and next; 
			$ta[$i] =~ m!^(\s*,*\s*)*$! and next; 
			my $tax_name = $hh[$i]; 
			$taxaN ++; 
			for my $tax_gene (split(/, /, $ta[$i])) {
				$geneN ++; 
				$back[-1]{'tax2gene'}{$tax_name}{$tax_gene} = $geneN; 
				defined $gene2grpID{$tax_name}{$tax_gene} and &stopErr("[Err] repeat ID {$tax_name}{$tax_gene}\n"); 
				$gene2grpID{$tax_name}{$tax_gene} = $back[-1]{'grpID'}; 
			}
		}
		$back[-1]{'geneN'} = $geneN; 
		$back[-1]{'taxaN'} = $taxaN; 
	}
	close($fh); 
	return(\@back, \%gene2grpID); 
}# load_grps_fmt_orthofinder() 

sub load_grps_fmt_orthocml {
	my $fn = shift; 
	my $fh = &openFH($fn, '<'); 
	my @back; 
	my %gene2grpID; 
	# For all_orthomcl.out 
	while (<$fh>) {
		chomp; 
		my @ta=split(/\t/, $_); 
		$ta[0] =~ m!^(ORTHOMCL\d+)\s*\(\s*(\d+)\s*genes\,\s*(\d+)\s*taxa\)\s*:\s*$! or die "$ta[0]\n"; 
		my $grpID = $1; 
		my $geneN = $2; 
		my $taxaN = $3; 
		push(@back, {}); 
		$back[-1]{'grpID'} = $grpID; 
		$back[-1]{'geneN'} = $geneN; 
		$back[-1]{'taxaN'} = $taxaN; 
		$ta[1] =~ s!^\s+|\s+$!!g; 
		my $i=0; 
		for my $tb (split(/\s+/, $ta[1])) {
			$i ++; 
			$tb =~ m!^(\S+)\(([^\s()]+)\)$! or die "|$tb|\n"; 
			my $tax_gene = $1; 
			my $tax_name = $2; 
			$back[-1]{'tax2gene'}{$tax_name}{$tax_gene} = $i; 
			defined $gene2grpID{$tax_name}{$tax_gene} and &stopErr("[Err] repeat ID {$tax_name}{$tax_gene}\n"); 
			$gene2grpID{$tax_name}{$tax_gene} = $back[-1]{'grpID'}; 
		}
	}
	close($fh); 
	return(\@back, \%gene2grpID); 
}# load_grps_fmt_orthocml () 

=head1 chk_prev_has( 'grpID|taxGene', [1,2,3,...], $grpID or [$taxID, $genID] )

Function : Check if gene or orthoID is in the previous categories; 

Return   : 1 - having ; 0 - not having

=cut
sub chk_prev_has {
	my ( $type, $prevN, $id ) = @_; 
	my @prev_k = map { "c${_}_has_$type" } @$prevN; 
	defined $id or &stopErr("[Err] id not defined\n"); 
	
	if ( $type eq 'grpID' ) {
		for my $k1 ( @prev_k ) {
			defined $glob{$k1} or &stopErr("[Err] Unknown k=[$k1]\n"); 
			defined $glob{$k1}{$id} and return (1); 
		}
	} elsif ( $type eq 'taxGene' ) {
		for my $k1 ( @prev_k ) {
			defined $glob{$k1} or &stopErr("[Err] Unknown k=[$k1]\n"); 
			defined $glob{$k1}{$id->[0]}{$id->[1]} and return (1); 
		}
	} else {
		&stopErr("[Err] Bad type [$type]\n"); 
	}
	return(0); 
}# sub chk_prev_has () 

=head1 load_taxaLis ($filename)

Return   : ( \%taxaLis )

  $taxaLis{ $taxaID } = $lineN 

=cut
sub load_taxaLis {
	my $fn = shift; 
	my %back; 
	my $fh = &openFH($fn, '<'); 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		$back{$ta[0]} //= $.; 
	}
	close($fh); 
	return(\%back); 
}# load_taxaLis () 

=head1 load_ggList () 

Function : Load all_12_taxa.gg 
Return   : (%gglis)
  $gglis{'tax2gene_h'} = { taxID=>{ geneID=>1 }, ... }
  $gglis{'tax2gene_a'} = { taxID=>[ geneID_1, geneID_2, ... ], ... }

=cut
sub load_ggList {
	my $fn = shift; 
	my %back; 
	my $fh = &openFH($fn, '<'); 
	while (&wantLineC($fh)) {
		my @ta = split(/\s+/, $_); 
		$ta[0] =~ m!^(\S+):$! or &stopErr("[Err] bad taxID [$ta[0]]\n"); 
		my $tax = $1; 
		for (my $i=1; $i<@ta; $i++) {
			$back{'tax2gene_h'}{$tax}{$ta[$i]} = $i; 
		}
		$back{'tax2gene_a'}{$tax} = [@ta[1..$#ta]]; 
		scalar( @{$back{'tax2gene_a'}{$tax}} ) == scalar( keys %{$back{'tax2gene_h'}{$tax}} ) or &stopErr("[Err] There might be repeat in [$ta[0]]\n"); 
	}
	close($fh); 
	return(\%back); 
}# load_ggList () 

=head1 out_data( $cN )


=cut
sub out_data {
	my ($cN) = @_; 
	( defined $glob{'odir'} and $glob{'odir'} ne '' ) or &stopErr("[Err] Bad -odir\n"); 
	$glob{'odir'} eq '.' and return; 
	my $k_grpID   = "c${cN}_has_grpID"; 
	my $k_taxGene = "c${cN}_has_taxGene"; 
	my %fh; 
	# grpID
	$fh{'grpID'} = &openFH("$glob{'odir'}/c${cN}.grpID", '>'); 
	for my $grpID (sort keys %{$glob{$k_grpID}}) {
		print {$fh{'grpID'}} "$grpID\n"; 
	}
	close( $fh{'grpID'} ); 

	# taxGene
	for my $tax ( sort keys %{$glob{$k_taxGene}} ) {
		$fh{'taxGene'} = &openFH("$glob{'odir'}/c${cN}.$tax.taxGene", '>'); 
		print {$fh{'taxGene'}} join("\t", qw/geneID taxID grpID/)."\n"; 
		for my $gene ( sort keys %{$glob{$k_taxGene}{$tax}} ) {
			my $grpID = 'NA'; 
			defined $taxGen2grpID{$tax}{$gene} and $grpID = $taxGen2grpID{$tax}{$gene}; 
			print {$fh{'taxGene'}} join("\t",$gene, $tax, $grpID )."\n"; 
		}
		close( $fh{'taxGene'} ); 
	}

	return; 
}# out_data

sub prepare_glob {
	for my $cN (qw/1 2 3 4 5 6 7 8 9/) {
		for my $p2 (qw/_has_taxGene _has_grpID _has_geneN/) {
			my $k = "c${cN}$p2"; 
			$glob{$k} = {}; 
		}
	}
	$glob{'odir'} = $opts{'odir'}; 
	if ($glob{'odir'} ne '.') {
		-d $glob{'odir'} and &stopErr("[Err] -odir $glob{'odir'} exists.\n"); 
		mkdir($glob{'odir'}) or &stopErr("[Err] Failed to create output_dir [$glob{'odir'}]"); 
	}
	$glob{'in_fmt'} = 'orthomcl'; 
	defined $opts{'in_fmt'} and $glob{'in_fmt'} = $opts{'in_fmt'}; 
}# prepare_glob

=head1 out_cN_geneN ([qw/1 2 3 4 5 6 7 8 9/])

Function : print gene numbers in each required $cN and taxa; 

Output format : 
category_ID \\t taxa_01 \\t taxa_02 \\t ... 
c1          \\t geneNum \\t geneNum \\t ...
c2          \\t geneNum \\t geneNum \\t ... 
...
All         \\t geneNum \\t geneNum \\t ...

=cut
sub out_cN_geneN {
	my ($cN_r) = @_; 
	my %taxID; 
	for my $k (keys %glob) {
		$k =~ m!^c(\d+)_has_geneN$! or next; 
		map { $taxID{$_} = 1; } keys %{$glob{$k}}; 
	}
	my @taxa_ids = sort keys %taxID; 
	print STDOUT join("\t", 'category_ID', @taxa_ids)."\n"; 
	for my $cN (@$cN_r) {
		my $k = "c${cN}_has_geneN"; 
		for my $tax (@taxa_ids) {
			$glob{$k}{$tax} //= 0; 
		}
		print STDOUT join("\t", "c$cN", @{$glob{$k}}{@taxa_ids})."\n"; 
	}
	my %ttl_taxGene_N; 
	for my $tax (@taxa_ids) {
		$ttl_taxGene_N{$tax} = 0; 
		defined $glob{'c4_gg_list'}{'tax2gene_a'}{$tax} and $ttl_taxGene_N{$tax} = scalar( @{$glob{'c4_gg_list'}{'tax2gene_a'}{$tax}} ); 
	}
	print STDOUT join("\t", "All", @ttl_taxGene_N{@taxa_ids})."\n"; 
}# out_cN_geneN ()


