package ReadInAlnSunhh; 
# Author : Honghe Sun, hs738@cornell.edu 
# Date   : 2014-07-18

use strict; 
use warnings; 
use LogInforSunhh; 
use Exporter qw(import);

our @EXPORT = qw(readMAF splitMafSline); 
our @EXPORT_OK = qw(normMAFloc); 

=head1 readMAF( $file_handle ) 

# MAF format : multiple alignment format. 
# https://cgwb.nci.nih.gov/FAQ/FAQformat.html#format5
Return       : (\%{ "a"=>[$Score_line1, $Score_line2, ... ], "o"=>[$Aligned_sequence_line1, $Aligned_sequence_line2, ... ] })
  'o' : Including ^'s '  and ^'p ' or something following ^'a ' 

Description  : Read in one record each time. 

Example of MAF format : 

 a score=9381 mismap=1e-09
 s Cmo_Scf00001  14516 2214 + 11258782 ACTGTAATTCTTCATTCTtttttttttttttttttggggggggggggggATCAAGAAATGGCATACAATGGATGAACTAAATACAATGCTACAGTAATTCAATATTCACAAAAGCGACT
 s Cma_Scf00047 710295 2197 +  1095042 ACTGTAATTCTTCATTCTttttttttttttttttttttgggggggggg---CAAGAAATGGTAGACAATGGATTTACTAAATACAATGCTACAGT----------TCACAAAAGCAACT
 p                                     "%+06<BHNTZ`flrx~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~xsmmmmmmmmmmsx~~~~~~~~~~~~
 p                                     &,28>DJPV\bhntz~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~xsmmmmmmmmmmrx~~~~~~~~~~~~

=cut
sub readMAF {
	my $fh = shift; 
	my %back_record; 
	my $is_read = 0; 
	while (<$fh>) {
		if ( m/^#/ ) {
		} elsif ( m/^$/ ) {
			# defined $back_record{a} or return {}; 
			defined $back_record{a} or next; 
			return \%back_record; 
		} elsif ( m/^a\s+/ ) {
			if ( $is_read == 0 ) {
				push(@{$back_record{a}}, $_); 
				$is_read = 1; 
			} else {
				&tsmsg("[Err] This is a bad format of MAF, because I cannot find an empty line to separate alignment blocks.\nLINE: $_\n"); 
				last; 
			}
		} else {
			if ( $is_read == 1 ) {
				push(@{$back_record{o}}, $_); 
			} else {
				&tsmsg("[Err] This is a bad format of MAF, because I get to \^s line before \^a line.\nLINE:$_\n"); 
			}
		}
	}
	defined $back_record{a} or return {}; 
	return \%back_record; 
}# sub readMAF()

=head1 splitMafSline( $MafSline_txt, $want_normal_position_[0|1] )

  $MafSline_txt         : Aligned "s" line from MAF format file. 
    Example : 
      s Cmo_Scf00001      0 2005 + 11258782 AATGTGAAACACCAACTGTCACACTGCACAGAACTGGATCCTGGAAACGGGGGAAGCAGAGTAGTGAGATAGCCTTTCCTTGAGGA
  $want_normal_position : Default 0. If 1, the two keys qw/normS normE/ were calculated as the real (1-based) position of start and end loci. 

Return       : (\%)
  \% : Keys = qw/seqId         seqStart blkSize seqStrand seqLen    seqSeq normS normE/
                 Cmo_Scf00001  0        2005    +         11258782  $AAT.  -1    -1
    (normS, normE) are all 1-based and normS < normE ; If $want_normal_position == 0, these two will be (-1, -1); 

=cut
sub splitMafSline {
	my $line = shift; 
	my $want_normP = shift; 
	( defined $want_normP and $want_normP ne '' ) or $want_normP = 0; 
	$line =~ s/\s+$//; 
	my @ta = split(/\s+/, $line); 
	$ta[0] eq 's' or die "Bad line: $line\n"; 
	my %back; 
	$back{seqId} = $ta[1]; 
	$back{seqStart} = $ta[2]; 
	$back{blkSize} = $ta[3]; 
	$back{seqStrand} = $ta[4]; 
	$back{seqLen} = $ta[5]; 
	$back{seqSeq} = $ta[6]; 
	$back{normS} = -1; 
	$back{normE} = -1; 
	if ($want_normP) {
		my @tse = &normMAFloc( $back{seqStrand}, $back{seqLen}, [$back{seqStart}, $back{blkSize}] ); 
		( $back{normS}, $back{normE} ) = @{$tse[0]}; 
	}
	return (\%back); 
}#End sub splitMafSline()

=head1 normMAFloc( $strand (+/-), $seqLen (number), @mafBlks )

  @mafBlks = ([mafStart_1, mafBlkSize_1], [mafStart_2, mafBlkSize_2], ...)

Return       : ( @normal_blkSE )
  @norm_blkSE = ([1-based-Block_Start, 1-based-Block_End], [1-based-Block_Start, 1-based-Block_End], ...)

=cut
sub normMAFloc {
	my ($str, $ttl_len, @mafBlks) = @_; 
	( $str eq '+' or $str eq '-' ) or &stopErr("[Err] Input Strand is [$str] instead of +/-\n"); 
	my @norm_blkSE; 
	if ( $str eq '+' ) {
		for (@mafBlks) {
			push(@norm_blkSE, [ $_->[0]+1, $_->[0]+$_->[1] ]); 
		}
	} else {
		for (@mafBlks) {
			my $blkE = $ttl_len - $_->[0]; 
			my $blkS = $blkE - $_->[1] + 1; 
			( $blkS > 0 and $blkE > 0 ) or &tsmsg("[Wrn]Converted normal (Start, End)=($blkS, $blkE) are not >= 1\n"); 
			push( @norm_blkSE, [$blkS, $blkE] ); 
		}
	}
	return (@norm_blkSE); 
}#End sub normMAFloc() 

1;
