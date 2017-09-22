#!/usr/bin/perl
### Convert maker gff3 to augustus format, and select gene models with one 5'-UTR and one 3'-UTR for the aug_UTR training. 
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
use fileSunhh; 

my %opts; 

GetOptions(\%opts, 
	"help!", 
	"utr3_n:i", "utr5_n:i", 
	"outFile:s", 
); 

sub usage {
	print <<UU;
##########################################################################################
# perl $0 in.maker_good.gff3
#
# -utr5_n         [N] Only genes with N 5'-UTR element will be output. (-1 as no checking)
# -utr3_n         [N] Only genes with N 3'-UTR element will be output. (-1 as no checking)
# 
# -outFile        [filename] file to create and write. 
##########################################################################################
UU
	exit 1; 
}

-t and !@ARGV and &usage(); 
$opts{'help'} and &usage(); 

my $oFH = ( defined $opts{'outFile'} ) ? &openFH($opts{'outFile'}, '>') : \*STDOUT ; 
$opts{'utr5_n'} = $opts{'utr5_n'} // -1; 
$opts{'utr3_n'} = $opts{'utr3_n'} // -1; 

my %fmt_name = qw(
CDS             CDS
THREE_PRIME_UTR 3'-UTR
3'-UTR          3'-UTR
FIVE_PRIME_UTR  5'-UTR
5'-UTR          5'-UTR
); 

my @geneLines; 
my $gID = ''; 
my %m2g ; 
while (<>) {
	m/^#/ and next; 
	m/^\s*$/ and next; 
	m/^>/ and last; 
	chomp; 
	my @ta = split(/\t/, $_); 
	if ( $ta[2] =~ m/^mRNA$/i or $ta[2] eq 'protein_match' ) {
		if ( scalar(@geneLines) > 0 ) {
			&writeGood($oFH, \@geneLines, $opts{'utr5_n'}, $opts{'utr3_n'}); 
			@geneLines = (); 
			$gID = ''; 
		}
		$ta[8] =~ m!^ID=([^\s;]+)! or &stopErr("[Err] $ta[8]\n"); 
		my $mID = $1; 
		$ta[8] =~ m!Parent=([^\s;]+)! or &stopErr("[Err] $ta[8]\n"); 
		$gID = $1; 
		defined $m2g{$mID} and &stopErr("[Err] Repeat mID=$mID\n"); 
		$m2g{$mID} = $gID; 
		$ta[8] = "transcript_id \"$mID\" ; gene_id \"$gID\";"; 
		# push(@geneLines, [[@ta], '']); 
	} elsif ( $ta[2] =~ m/^(?:cds|five_prime_UTR|three_prime_UTR)$/i or $ta[2] eq 'match_part' ) {
		$ta[8] =~ m!Parent=([^\s;]+)! or &stopErr("[Err] No parent : $_\n"); 
		my $pIDs = $1; 
		my @transIDs; 
		my $tgID; 
		for my $pID (split(/,/, $pIDs)) {
			push(@transIDs, "transcript_id \"$pID\";"); 
			defined $m2g{$pID} or &stopErr("[Err] Unkown m2g{$pID}\n"); 
			$tgID = $tgID // $m2g{$pID}; 
			$tgID eq $m2g{$pID} or &stopErr("[Err] parent different: $_\n"); 
		}
		$ta[8] = join(" ", @transIDs, "gene_id \"$tgID\""); 
		
		$ta[2] = uc($ta[2]); 
		defined $fmt_name{$ta[2]} and $ta[2] = $fmt_name{$ta[2]}; 

		defined $fmt_name{$ta[2]} and push(@geneLines, [[@ta], '']); 
		
	} elsif ( $ta[2] =~ m/^gene$/i ) {
		$gID = '' ; 
		$ta[8] =~ m/^ID=([^\s;]+)/ or &stopErr("[Err] Failed to parse gene element : $_\n"); 
		$gID = $1; 
	} else {
		$ta[2] =~ m/^exon$/i and next; 
		&stopErr("[Err] 1: $_\n"); 
	}
}

if ( scalar(@geneLines) > 0 ) {
	&writeGood($oFH, \@geneLines, $opts{'utr5_n'}, $opts{'utr3_n'}); 
	@geneLines = (); 
	$gID = ''; 
}

sub writeGood {
	my ($fh, $ar, $n_utr5, $n_utr3) = @_; 
	$n_utr5 = $n_utr5 // -1; 
	$n_utr3 = $n_utr3 // -1; 

	my $is_good = 1; 
	
	my ($n5, $n3)=(0,0); 
	for (my $i=0; $i<@$ar; $i++) {
		$ar->[$i][0][8] =~ m/transcript_id \".*transcript_id \"/ and do { $is_good = 0; last; }; 
		$n_utr5 != -1 and $ar->[$i][0][2] eq "5'-UTR" and $n5++; 
		$n_utr3 != -1 and $ar->[$i][0][2] eq "3'-UTR" and $n3++; 
	}
	$n_utr5 != -1 and $n5 != $n_utr5 and $is_good = 0; 
	$n_utr3 != -1 and $n3 != $n_utr3 and $is_good = 0; 
	
	$is_good == 0 and return ; 
	
	for (my $i=0; $i<@$ar; $i++) {
		print {$fh} join("\t", @{$ar->[$i][0]})."\n";
	}
	return 0;
}# sub writeGood() 

