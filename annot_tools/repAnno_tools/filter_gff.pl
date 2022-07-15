#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

# ###
# seq999  LTRharvest      repeat_region   195972  198652  .       +       .       ID=repeat_region813
# seq999  LTRharvest      target_site_duplication 195972  195976  .       +       .       Parent=repeat_region813
# seq999  LTRharvest      inverted_repeat 195977  195978  .       +       .       Parent=repeat_region813
# seq999  LTRharvest      LTR_retrotransposon     195977  198647  .       +       .       ID=LTR_retrotransposon813;Parent=repeat_region813;ltr_similarity=99.70;seq_numbe
# seq999  LTRharvest      long_terminal_repeat    195977  196312  .       +       .       Parent=LTR_retrotransposon813
# seq999  LTRdigest       primer_binding_site     196315  196325  5.44    +       .       Parent=LTR_retrotransposon813;trna=Cryptococcus_neoformans_var_JEC21_chr4.trna5-
# seq999  LTRharvest      long_terminal_repeat    198312  198647  .       +       .       Parent=LTR_retrotransposon813
# seq999  LTRharvest      inverted_repeat 196311  196312  .       +       .       Parent=repeat_region813
# seq999  LTRharvest      inverted_repeat 198312  198313  .       +       .       Parent=repeat_region813
# seq999  LTRharvest      inverted_repeat 198646  198647  .       +       .       Parent=repeat_region813
# seq999  LTRharvest      target_site_duplication 198648  198652  .       +       .       Parent=repeat_region813

my $max_dist_to_ltr = 20; 
my $minR_inInclRegion = 0.5; 

print STDOUT join("\t", qw/seqID LTR1_S LTR1_E LTR2_S LTR2_E Inner_S Inner_E PBS_S PBS_E Strand/)."\n"; 
my %rec; 
my @id_list; 
while (<>) {
	m/^##s/ and next; 
	m/^#S/ and next; 
	m/^#/ and next; 
	chomp; 
	my @ta = split(/\t/, $_); 
	if ( $ta[1] eq 'LTRharvest' and $ta[2] eq 'repeat_region' ) {
		%rec = (); 
	} elsif ( $ta[1] eq 'LTRharvest' and $ta[2] eq 'long_terminal_repeat' ) {
		my $la2 = <>; chomp($la2); 
		my @ta2 = split(/\t/, $la2); 
		if ( $ta2[1] eq 'LTRdigest' ) {
			$ta2[2] eq 'primer_binding_site' or die "$_\n"; 
			my $la3 = <>; chomp($la3); 
			my @ta3 = split(/\t/, $la3); 
			( $ta3[1] eq 'LTRharvest' and $ta3[2] eq 'long_terminal_repeat' ) or die "$_\n"; 
			my ($ltr_s, $ltr_e) = ( $ta[3], $ta3[4] ); 
			my ($inner_s, $inner_e) = ( $ta[4]+1, $ta3[3]-1 ); 
			my ($pbs_s, $pbs_e) = ($ta2[3], $ta2[4]); 
			&olapLen( $inner_s, $inner_e, $pbs_s, $pbs_e ) >= $minR_inInclRegion * ($pbs_e-$pbs_s+1) or next; 
			&minDist( $inner_s, $inner_e, $pbs_s, $pbs_e ) <= $max_dist_to_ltr or next; 
			print STDOUT join("\t", $ta[0], $ta[3], $ta[4], $ta3[3], $ta3[4], $inner_s, $inner_e, $pbs_s, $pbs_e, $ta[6])."\n"; 

		} elsif ( $ta2[1] eq 'LTRharvest' ) {
			$ta2[2] eq 'long_terminal_repeat' or die "$_\n"; 
		} else {
			die "$_\n"; 
		}
	}
}

sub olapLen {
	my ($s1,$e1, $s2, $e2) = @_; 
	my $s3 = ($s1 > $s2) ? $s1 : $s2 ; 
	my $e3 = ($e1 < $e2) ? $e1 : $e2 ; 
	return ($e3-$s3+1) ; 
}

sub minDist {
	my ($s1,$e1, $s2, $e2) = @_; 
	$e2 < $s1 and die "$e2 < $s1 : ($s1,$e1, $s2, $e2)\n"; 
	$s2 > $e1 and die "$s2 > $e1 : ($s1,$e1, $s2, $e2)\n"; 
	my $dist1 = $s2-$s1-1; 
	my $dist2 = $e2-$e1-1; 
	my $dist = ($dist1 > $dist2) ? $dist2  : $dist1 ; 
	return $dist ; 
}



