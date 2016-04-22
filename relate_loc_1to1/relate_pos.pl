#!/usr/bin/perl 
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use ReadInAlnSunhh; 

-t and !@ARGV and die "perl $0 in_1to1.maf > in_1to1.maf.rel_loc\n"; 

my $inFh = \*STDIN; 
@ARGV and $inFh = &openFH( shift, '<' ); 

my %cnt; 
$cnt{'cntN_step'} = 5e4; 
while ( my $rec = &ReadInAlnSunhh::readMAF($inFh) ) {
	&fileSunhh::log_section( $. , \%cnt ) and &tsmsg("[Msg] Processing $. line.\n"); 
	(keys %$rec) > 0 or last; 
	defined $rec->{'o'} or next; 
	my $h1 = &ReadInAlnSunhh::splitMafSline( $rec->{'o'}[0], 1 ); 
	my $h2 = &ReadInAlnSunhh::splitMafSline( $rec->{'o'}[1], 1 ); 
	my @seq_a1 = split(//, $h1->{'seqSeq'}); 
	my @seq_a2 = split(//, $h2->{'seqSeq'}); 

	my ($p1, $p2) = (0, 0); 
	my @p2p; 
	for (my $i=0; $i<@seq_a1; $i++) {
		if ( $seq_a1[$i] =~ m/^[ATGCN]$/i ) {
			if ( $seq_a2[$i] =~ m/^[ATGCN]$/i ) {
				$p1 ++; $p2 ++; 
			} else {
				$p1 ++; 
			}
		} else {
			if ( $seq_a2[$i] =~ m/^[ATGCN]$/i ) {
				$p2 ++; 
			} else {
				&stopErr("[Err] Bad input at [$p1, $p2] with [$seq_a1[$i], $seq_a2[$i]]\n"); 
			}
		}

		push(@p2p, [$p1, $p2]); 
	}

	if ( $h1->{'seqStrand'} eq '+' ) {
		if ( $h2->{'seqStrand'} eq '+' ) {
			@p2p = map { [$_->[0]+$h1->{'normS'}-1, $_->[1]+$h2->{'normS'}-1, '+'] } @p2p; 
		} else {
			@p2p = map { [$h1->{'normE'}-$_->[0]+1, $h2->{'normE'}-$_->[1]+1, '-'] } @p2p; 
		}
	} else {
		&stopErr("[Err] I don't like this format! Strand of first seq is not '+': $rec->{'o'}[0]\n"); 
	}

	for my $tp (@p2p) {
		print STDOUT join("\t", $h1->{'seqId'}, $tp->[0], $h2->{'seqId'}, $tp->[1], $tp->[2])."\n"; 
	}

	print STDOUT "\n"; 
}

&tsmsg("[Rec] All done for $0\n"); 
# a score=9381 mismap=1e-09
# s Cmo_Scf00001  14516 2214 + 11258782 ACTGTAATTCTTCATTCTtttttttttttttttttggggggggggggggATCAAGAAATGGCATACAATGGATGAACTAAATACAATGCTACAGTAATTCAATATTCACAAAAGC
# s Cma_Scf00047 710295 2197 +  1095042 ACTGTAATTCTTCATTCTttttttttttttttttttttgggggggggg---CAAGAAATGGTAGACAATGGATTTACTAAATACAATGCTACAGT----------TCACAAAAGC
# p                                     "%+06<BHNTZ`flrx~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~xsmmmmmmmmmmsx~~~~~~~~

