#!/usr/bin/perl
# 2014-07-16 Design to use multiple flank lengths to overcome disadvantage of overlap alignment method. 
use strict; 
use warnings; 
use SeqAlnSunhh; # For olap_e2e_A2B( seqA, seqB [, {para_hash}])
use Getopt::Long; 

my %opts; 
GetOptions(\%opts, 
	'help!', 
	"maxNlen:i", # length of max N gap. 
	"flankLen:s", # length of maximum flank region to be checked for overlapping 
	"minIdent:f", # Minimum identity to accept an overlap alignment
	"minAlnLen:i", # Minimum length of overlap alignment required. 
	"opt_tail:i", # Allow trimming $opts{opt_tail} bp from tail of the upstream subseq (seqA) in the overlap alignment. 
	"opt_head:i", # Allow trimming $opts{opt_head} bp from head of the downstream subseq (seqB) in the overlap alignment. 
	"logChange:s", 
	"outSeqF:s", 
	"maxSingleRatio:f", 
); 

-t and !@ARGV and &usage(); 
$opts{hep} and &usage(); 

defined $opts{maxNlen} or $opts{maxNlen} = 1; 
defined $opts{flankLen} or $opts{flankLen} = 200; 
defined $opts{minIdent} or $opts{minIdent} = 0.99; 
defined $opts{minAlnLen} or $opts{minAlnLen} = 15; 
defined $opts{maxSingleRatio} or $opts{maxSingleRatio} = 0.9;
my $para_olap = {}; 
defined $opts{opt_tail} and $para_olap->{opt_tail} = $opts{opt_tail}; 
defined $opts{opt_head} and $para_olap->{opt_head} = $opts{opt_head}; 

sub usage {
	print STDOUT <<HELP; 
################################################################################
# perl $0 to_chk.scf.fa 
# -help
# 
# -outSeqF     File to store changed sequences. [*STDOUT]
# -logChange   File to record changed regions. [*STDERR]
#
# -maxNlen     Max length of N gap to check. [1]
# -flankLen    Max lengths of flank region to be used for overlapping alignment[200[,50,100,200]].
# -minIdent    Min identity to accept an overlap alignment[0.99]
# -minAlnLen   Min length to accept an overlap alignment[15]
# -opt_tail    Max Bp to remove from tail of upstream sub-sequence for alignment.[0]
# -opt_head    Max Bp to remove from head of downstream subSeq for alignment.[0]
#
################################################################################
HELP
	exit 1; 
}

my @FlankLen; 
for (split(/,/, $opts{flankLen})) {
	s/\s//g; 
	push(@FlankLen, $_); 
}
my $logChangeFH = \*STDERR; 
if (defined $opts{logChange}) {
	undef($logChangeFH); 
	open ($logChangeFH, ">", "$opts{logChange}") or die ; 
}
print {$logChangeFH} join("\t", qw/ScafID RawStart RawEnd RemovedSeq Chk:Seg1:Seg2 Aln:Blk1:Blk2 Aln_size Aln_ident Aln_score Aln_Seqs/)."\n"; 

my $outSeqFH = \*STDOUT; 
if (defined $opts{outSeqF}) {
	undef($outSeqFH); 
	open ($outSeqFH, ">", "$opts{outSeqF}") or die ; 
}


## Read in sequences; 
my (%seq, @IDs); 
{
	&tsmsg("[Rec] Begin to read in sequences.\n"); 
	my $kk;
	while (<>) {
		if ( m/^\s*>(\S+)/ ) {
			$kk = $1;
		} else {
			defined $seq{$kk} or push(@IDs, $kk);
			$seq{$kk} .= $_;
		}
	}
	my $seq_num = $#IDs+1; 
	&tsmsg("[Rec] Total [$seq_num] sequences to check.\n"); 
}

my $seqN = 0; 
for my $tk (@IDs) {
	my @change_info; 
	# Load scaffold sequence. 
	$seqN ++; 
	$seq{$tk} =~ s/\s//g; 
	my $ll = length($seq{$tk}); 
	&tsmsg("[Msg] Dealing with [$seqN] $tk len=$ll\n"); 

	# Fill N gaps along scaffold one by one. 
	my $newseq = ''; 
	pos($seq{$tk}) = 0; 
	my $prevE = -1; # Position in $seq{$tk} for the last base from $newseq 
	MMM: 
	while ( $seq{$tk} =~ m/\G(?:.*?)([nN]+)/gs ) {
		my $n_seq = $1; 
		my ($gapS, $gapE) = ($-[1]+1, $+[1]); 
		pos($seq{$tk}) = $gapE; 
		if ($newseq eq '') {
			if ( $gapS > 1 ) {
				$newseq = substr($seq{$tk}, 0, $gapS-1); 
				$prevE = $gapS-1; 
			} else {
				# The N is in the head of sequence. I do not want to remove it. So just record the N gap sequence. 
				$newseq = substr($seq{$tk}, 0, $gapE); 
				$prevE = $gapE; 
				next MMM; 
			}
		} else {
			if ( ($gapS-1) - ($prevE+1) + 1 > 0 ) {
				$newseq = $newseq . substr( $seq{$tk}, $prevE, ($gapS-1) - ($prevE+1) + 1 ); 
				$prevE = $gapS-1; 
			} else {
				# $newseq does not need to be changed. 
				; 
			}
		}
		my $is_removed = 0; 
		if ( $gapE-$gapS+1 <= $opts{maxNlen} ) {
			CHK_LEN: 
			for my $flankLen (@FlankLen) {
				# Extract subseq from $newseq. 
				my ($seqA) = ( $newseq =~ m/([ATGC]{1,$flankLen})$/i ); 
				my $lenA = length($seqA); 
				if ( $lenA >= $opts{minAlnLen}) {
					my $old_sA = $prevE-$lenA+1; # This is not quite exactly right, especially when there are two close N gaps. 
					my $old_eA = $prevE; 
					my $seqB = substr( $seq{$tk}, $gapE, $flankLen ); 
					$seqB =~ s/[nN].*$//; 
					my $lenB = length($seqB); 
					if ( $lenB >= $opts{minAlnLen} ) {
						my $old_sB = $gapE+1; 
						my $old_eB = $gapE+1+$lenB-1; 
						# Now we should check overlap in this case. 
						my %aln = %{ &olap_e2e_A2B($seqA, $seqB, $para_olap) }; 
						if ( $aln{aln_ident} >= $opts{minIdent} and $aln{endA}-$aln{startA}+1 >= $opts{minAlnLen} and $aln{endB}-$aln{startB}+1 >= $opts{minAlnLen}) {
							&mainBaseRatio($aln{seqA_aln}) >= $opts{maxSingleRatio} and next CHK_LEN; 
							# Here we choose to use the upstream sequences. 
							my $rmTailAlen = $lenA-$aln{endA}; 
							if ( $rmTailAlen > 0 ) {
									substr($newseq, -$rmTailAlen) = ''; 
							}
							$prevE = $old_sB + $aln{endB}-1; # Update end position on old sequence. $prevE
							# Record positioins. 
							push(@change_info, [
								$tk, # Old sequence name. 
								$old_eA-$rmTailAlen+1, # Start position to remove from old sequence. 
								$prevE, # End position to remove from old sequence. 
								substr($seqA, $aln{endA}) . $n_seq . substr($seqB, 0, $aln{endB}) , # Sequence to remove from old sequence. 
								join(':', "Chk", "($old_sA-$old_eA)[$lenA]", "($old_sB-$old_eB)[$lenB]"), 
								join(':', "Aln", "($aln{startA}-$aln{endA})", "($aln{startB}-$aln{endB})"), 
								$aln{aln_len}, 
								$aln{aln_ident}, 
								$aln{aln_score}, 
								join(":", $aln{seqA_aln}, $aln{seqB_aln}, $aln{seqC_aln})
							]); 
							$is_removed = 1; 
							last CHK_LEN; 
						}
					} else { # lenB is too small to do alignment. 
					}
				} else { # lenA is too small to do alignment. 
				}
			}#End for CHK_LEN; 
		} else { # This gap is too large, so do not need to check overlap. 
		}
		$is_removed == 0 and do { $newseq = $newseq . $n_seq; $prevE = $gapE; }; 
		pos($seq{$tk}) = $gapE; 
	}
	# Add the rest sequence if needed. 
	if ( $prevE == -1 ) {
		$newseq = $seq{$tk}; 
	} elsif ( $prevE < $ll ) {
		$newseq = $newseq . substr($seq{$tk}, $prevE); 
		$prevE = $ll; 
	}

	# Output merged sequence. $newseq
	&tsmsg("[Msg] Output new sequence of [$tk]\n"); 
	# my $new_ll = length($newseq); 
	print {$outSeqFH} ">$tk\n$newseq\n"; 
	for my $tr (@change_info) {
		print {$logChangeFH} join("\t", @$tr)."\n"; 
	}
}

&tsmsg("[Rec]All done.\n"); 


sub mainBaseRatio {
	my %cc; 
	$_[0] = uc($_[0]); 
	my $nn = 0; 
	while ( $_[0] =~ m/(.)/g ) {
		$cc{$1} ++; 
		$nn ++; 
	}
	my $vv = (sort { $b <=> $a } values %cc)[0]; 
	return $vv/$nn ; 
}

sub tsmsg {
	my $tt = scalar(localtime()); 
	print STDERR join('', "[$tt]", @_); 
}

