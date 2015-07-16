#!/usr/bin/perl 
# Clip scaffolds' end to make them fit the following scaffolding alignment. 

use strict;
use warnings; 
use ReadInSeqSunhh; 
# use DealFastaSunhh; 
use fastaSunhh; 
use LogInforSunhh; 
use Getopt::Long; 

!@ARGV and die "perl $0 in.scf.fa\n"; 

my %opts; 
GetOptions(\%opts, 
	"chk_len:i", 
	"max_gap:i", 
	"min_noN:i", 
	"min_clust_len:i", 
	"min_outLen:i", 
	"help!"
); 

my $chk_len = (defined $opts{chk_len}) ? $opts{chk_len} : 10e3; 
my $max_gap = (defined $opts{max_gap}) ? $opts{max_gap} : 700; 
my $min_noN = (defined $opts{min_noN}) ? $opts{min_noN} : 100; 
my $min_clust_len = (defined $opts{min_clust_len}) ? $opts{min_clust_len} :  3e3; 
my $min_outLen = (defined $opts{min_outLen}) ? $opts{min_outLen} : 100; 

my $srh_pat = '[ATGCatgc]+'; 

my $seqF = shift; 
my $sfh; 
open $sfh,'<',"$seqF" or die; 
my $sc_ct = 0; 
for ( my ($relHR, $get) = &get_fasta_seq($sfh); defined $relHR and $get; ($relHR, $get) = &get_fasta_seq($sfh) ) {
	my $tmpSeq = $relHR->{seq}; 
	$tmpSeq =~ s/\s//g; 
	# From 5' 
	my ($has_good5, @clipped_5) = &clip_5p($tmpSeq, $chk_len, $max_gap, $min_noN, $min_clust_len); 
	my ($has_good3, @clipped_3); 
	if ( $has_good5 ) {
		my $restSeq = pop(@clipped_5); 
		$restSeq =~ tr/ATGCatgc/TACGtacg/; 
		$restSeq = reverse($restSeq); 
		($has_good3, @clipped_3) = &clip_5p($restSeq, $chk_len, $max_gap, $min_noN, $min_clust_len); 
		for my $ts (@clipped_3) {
			$ts =~ tr/ATGCatgc/TACGtacg/; 
			$ts = reverse($ts); 
		}
		push(@clipped_5, @clipped_3); 
	}
	if (scalar(@clipped_5) == 1) {
		# raw sequence not changed. 
		unless ( length ( $clipped_5[0] ) == length($tmpSeq) and $clipped_5[0] eq $tmpSeq ) {
			&tsmsg( "[Err] key=$relHR->{key} problem: Maybe terminated with Ns.\nclipped: $clipped_5[0]\nrawSeq : $tmpSeq\n" ); 
		}
		if ( length($clipped_5[0]) >= $min_outLen ) {
			my $ts = $clipped_5[0]; $ts =~ s/(.{100})/$1\n/g; chomp($ts); 
			print STDOUT ">$relHR->{head}\n$ts\n"; 
		}
	} else {
		for my $ts (@clipped_5) {
			$ts =~ s/^[nN]+//; $ts =~ s/[nN]+$//; 
			length($ts) >= $min_outLen or next; 
			$sc_ct ++; 
			$ts =~ s/(.{100})/$1\n/g; chomp($ts); 
			print STDOUT ">Clip.$sc_ct [From $relHR->{key}]\n$ts\n"; 
		}
	}
}
close $sfh; 

# Return : ($has_good, @clipped_seq)
sub clip_5p {
	my ($tmpSeq, $chk_len, $max_gap, $min_noN, $min_clust_len) = @_; 
	my @clipped_seq; 
	# From the 5' direction. 
	my $need_clip5 = 1; 
	my ($prevS, $prevE) = (-1, -1); # (start, end) of good block. 
	
	if ( length($tmpSeq) <= $chk_len ) {
		# Split sequence according to $max_gap
		push( @clipped_seq, split(/[nN]{$max_gap,}/, $tmpSeq) ); 
		$tmpSeq = ''; 
	} else {
		# Search good noN regions loci. 
		my @posList = &siteList(\$srh_pat, \$tmpSeq, 'Min'); 
		## Remove match (noN) region shorter than $min_noN ; 
		{
		my @newArr; 
		for my $ar0 (@posList) {
			$ar0->[1] - $ar0->[0] + 1 >= $min_noN and push(@newArr, $ar0); 
		}
		@posList = @newArr; 
		}
		
		# Search for good block. 
		BLK: 
		for ( my $i=0; $i<@posList and $need_clip5 == 1; $i++ ) {
			if ( $prevS == -1 ) {
				$posList[$i][0] > 1 and push( @clipped_seq, split(/[nN]{$max_gap,}/, substr($tmpSeq, 0, $posList[$i][0]-1)) ); 
				($prevS, $prevE) = ( $posList[$i][0], $posList[$i][1] ); 
				if ( $prevE-$prevS+1 >= $min_clust_len ) {
					$tmpSeq = substr($tmpSeq, $prevS-1); 
					$need_clip5 = 0; 
					last BLK; 
				}
			} elsif ( $posList[$i][0]-$prevE-1 <= $max_gap ) {
				# Should be joined here. 
				$prevE = $posList[$i][1]; 
				if ( $prevE-$prevS+1 >= $min_clust_len ) {
					$tmpSeq = substr($tmpSeq, $prevS-1); 
					$need_clip5 = 0; 
					last BLK; 
				}
			} elsif ( $posList[$i][0]-$prevE-1 > $max_gap ) {
				# Should be splitted here. And the previous block [$prevS, $prevE] is not long enough. 
				push( @clipped_seq, substr($tmpSeq, $prevS-1, $prevE-$prevS+1) ); 
				push( @clipped_seq, substr($tmpSeq, $prevE, $posList[$i][0]-$prevE+1) ); 
				( $prevS, $prevE ) = ( $posList[$i][0], $posList[$i][1] ); 
				if ( $prevE-$prevS+1 >= $min_clust_len ) {
					$tmpSeq = substr($tmpSeq, $prevS-1); 
					$need_clip5 = 0; 
					last BLK; 
				}
			} else {
				&stopErr("[Err] Why we are here!\n"); 
			}
		}#End for 
	}
	
	
	if ( $prevS == -1 ) { 
		# There isn't any good block in the input sequence. 
		$tmpSeq ne '' and push( @clipped_seq, split(/[nN]{$max_gap,}/, $tmpSeq) ); 
	} elsif ($need_clip5 == 0) {
		# The good block starts from the 5' of $tmpSeq. 
		push(@clipped_seq, $tmpSeq); 
	} else {
		# The good long enough block does not exist. 
		push( @clipped_seq, substr($tmpSeq, $prevS-1, $prevE-$prevS+1) ); 
		$tmpSeq = substr($tmpSeq, $prevE); 
		push( @clipped_seq, split(/[nN]{$max_gap,}/, $tmpSeq) ); 
	}
	
	my $has_good = ( $need_clip5 ) ? 0 : 1 ; 
	return($has_good, @clipped_seq); 
}#End sub clip_5p

