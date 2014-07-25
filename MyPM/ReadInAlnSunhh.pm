package ReadInAlnSunhh; 
# Author : Honghe Sun, hs738@cornell.edu 
# Date   : 2014-07-18

use strict; 
use warnings; 
use Exporter qw(import);

our @EXPORT = qw(readMAF splitMafSline); 
our @EXPORT_OK = qw(normMAFloc); 


# Input : (FILE_Handle)
# Return: \%{ "a"=>Score_lines, "o"=>Aligned_sequence_lines }
# MAF format : multiple alignment format. 
# https://cgwb.nci.nih.gov/FAQ/FAQformat.html#format5
sub readMAF {
	my $fh = shift; 
	my $curpos = tell($fh); 
	my %back_record; 
	my $is_read = 0; 
	while (<$fh>) {
		if ( m/^#/ ) {
		} elsif ( m/^a\s+/ ) {
			if ( $is_read == 0 ) {
				push(@{$back_record{a}}, $_); 
				$is_read = 1; 
			} else {
				seek( $fh, $curpos, 0 ); 
				last; 
			}
		} else {
			push(@{$back_record{o}}, $_); 
		}
		$curpos = tell($fh); 
	}
	defined $back_record{a} or return {}; 
	return \%back_record; 
}# sub readMAF()

# Input  : ($MafSline[, 0/1])
# Input  : Single "s" line from MAF format file. 
# Return : \%{} : Keys = qw/seqId seqStart blkSize seqStrand seqLen seqSeq normS normE/
#          (normS, normE) are all 1-based and normS < normE ; 
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
	return \%back; 
}#End sub splitMafSline()

# Input : ($strand, $total_length, @mafBlks)
#          @mafBlks = ([mafStart, mafBlkSize], [mafStart, mafBlkSize], ...)
# Output: (@normal_blkSE)
#          @norm_blkSE = ([1-based-Block_Start, 1-based-Block_End], [1-based-Block_Start, 1-based-Block_End], ...)
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

sub stopErr {
	&tsmsg(@_); 
	exit 1; 
}
sub tsmsg {
	my $tt = scalar(localtime()); 
	print STDERR join('', @_); 
}

1;
