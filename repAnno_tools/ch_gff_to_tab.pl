#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

-t and !@ARGV and die "perl $0 PG1All_v2.scf.fa.gff99.dgt 1> dgt.tab\n"; 
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

print STDOUT join("\t", qw/eleID eleS eleE Str seqID LTR1_S LTR1_E LTR2_S LTR2_E Inner_S Inner_E PBS_S PBS_E PPT_S PPT_E/)."\n"; 
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
		$ta[8] =~ m/^ID=([^\s;]+)$/ or &stopErr("parse $ta[8]\n"); 
		$rec{eleID} = $1; 
		$rec{eleID} =~ s/^repeat_region/RR/; 
		$rec{eleS} = $ta[3]; 
		$rec{eleE} = $ta[4]; 
		$rec{eleStrand} = $ta[6]; 
	} elsif ( $ta[1] eq 'LTRharvest' and $ta[2] eq 'long_terminal_repeat' ) {
		my %infor; 
		$infor{seqID} = $ta[0]; 
		$infor{ltr1S} = $ta[3]; 
		$infor{ltr1E} = $ta[4]; 
		$infor{innS}  = $infor{ltr1E}+1; 
		
		SUB_LINE: 
		while ( my $la2 = <> ) {
			chomp($la2); 
			my @ta2 = split(/\t/, $la2); 
			if ( $ta2[1] eq 'LTRdigest' ) {
				if ( $ta2[2] eq 'primer_binding_site' ) {
					$infor{pbsS} = $ta2[3]; 
					$infor{pbsE} = $ta2[4]; 
				} elsif ( $ta2[2] eq 'RR_tract' ) {
					$infor{pptS} = $ta2[3]; 
					$infor{pptE} = $ta2[4]; 
				} else {
					die "1:$la2\n"; 
				}
			} elsif ( $ta2[1] eq 'LTRharvest' and $ta2[2] eq 'long_terminal_repeat' ) {
				$infor{ltr2S} = $ta2[3]; 
				$infor{ltr2E} = $ta2[4]; 
				$infor{innE}  = $infor{ltr2S}-1; 
				last SUB_LINE; 
			} else {
				die "2:$la2\n"; 
			}
		}
		for my $t (qw/pbsS pbsE pptS pptE/) {
			defined $infor{$t} or $infor{$t} = -1; 
		}
		print STDOUT join("\t", @rec{qw/eleID eleS eleE eleStrand/},  @infor{qw/seqID ltr1S ltr1E ltr2S ltr2E innS innE pbsS pbsE pptS pptE/})."\n"; 
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



