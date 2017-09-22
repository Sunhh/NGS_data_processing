#!/usr/bin/perl
use strict; 
use warnings; 
use Cwd 'abs_path'; 
use LogInforSunhh; 
use fileSunhh; 
use gffSunhh; 
use mathSunhh; 
my $ms = mathSunhh->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"in_gff3:s", # input gff3 file 
	"out_gff3:s", # Output gff3 file 
	"scafKLfile:s", # [filename] key\tlength for scaffold.
	
	# Filter options 
	"trimOverflow!", # Remove scaffold regions overflowing the maximum length of scaffold.
	"rm_bad1st!",    # Remove spaln_prot2genome alignments' first bad alignment. 
	
	"printCmd!", 
	"help!", 
); 

!@ARGV and die "perl $0 CM_CS_WM97_Arab_SprotPln.spaln_M4.fmt.gff3\n"; 

my $inFh = \*STDIN; 
my $outFh = \*STDOUT; 
defined $opts{'in_gff3'} and $inFh = &openFH( $opts{'in_gff3'}, 'read' ); 
if ( defined $opts{'out_gff3'} ) {
	$outFh = &openFH( $opts{'out_gff3'}, 'write' ); 
}

my $gs = gffSunhh->new(); 
my ($gff_href) = $gs->read_gff3File(
 'gffFH'=>$inFh, 
 'top_hier'=>{ 'mrna'=>1, 'match'=>2, 'protein_match'=>3, 'expressed_sequence_match'=>4 }
); 
my %gff_hash = %$gff_href; 


my @topIDs = sort { $gff_hash{'lineN_group'}{$a}{'curLn'}[0] <=> $gff_hash{'lineN_group'}{$a}{'curLn'}[0] } keys %{$gff_hash{'lineN_group'}}; 
print {$outFh} join("\t", qw/TopID TopS TopE TopLen ExonNum AvgExonLen AvgIntronLen/)."\n"; 
for my $topID ( @topIDs ) {
	my ($topS, $topE, $topLen, $exonNum, $avgExonLen, $avgIntronLen); 
	($topS, $topE, $topLen) = (-1, -1, -1); 
	my $curLn = $gff_hash{'lineN_group'}{$topID}{'curLn'}[0]; 
	# Output parent lines. 
	for my $parLn (@{$gff_hash{'lineN_group'}{$topID}{'parLn'}}) {
		# defined $gff_hash{'lineN2line'}{$parLn} and print {$outFh} "$gff_hash{'lineN2line'}{$parLn}\n"; 
	}
	# Output current line. 
	# defined $gff_hash{'lineN2line'}{$curLn} and print {$outFh} "$gff_hash{'lineN2line'}{$curLn}\n"; 
	if ( defined $gff_hash{'lineN2line'}{$curLn} ) {
		$topS = $gff_hash{'lineN2hash'}{$curLn}{'start'}; 
		$topE = $gff_hash{'lineN2hash'}{$curLn}{'end'}; 
		$topLen = $topE-$topS+1; 
	}
	# Output offspring lines. 
	my %len; 
	my (@len_exon, @len_intron); 
	$len{'exon'} = []; 
	$len{'intron'} = []; 
	for my $offLn (@{$gff_hash{'lineN_group'}{$topID}{'offLn'}}) {
		defined $gff_hash{'lineN2line'}{$offLn} or next; 
		if ($gff_hash{'lineN2hash'}{$offLn}{'type'} =~ m/^match_part$/i) {
			my %t_hash = %{ $gff_hash{'lineN2hash'}{$offLn} }; 
			if ( @{$len{'exon'}} > 0 ) {
				if ($t_hash{'strand'} eq '+') {
					push(@{$len{'intron'}}, [ $t_hash{'start'}-$len{'exon'}[-1][2]-1 ]); 
					push(@len_intron, $len{'intron'}[-1][0]); 
				} elsif ($t_hash{'strand'} eq '-') {
					push(@{$len{'intron'}}, [ $len{'exon'}[-1][1]-$t_hash{'end'}-1 ]); 
					push(@len_intron, $len{'intron'}[-1][0]); 
				} else {
					die "$gff_hash{'lineN2line'}{$offLn}\n"; 
				}
			}
			push(@{$len{'exon'}}, [$t_hash{'end'}-$t_hash{'start'}+1, $t_hash{'start'}, $t_hash{'end'}]); 
			push(@len_exon, $len{'exon'}[-1][0]); 
		} else {
			die "$gff_hash{'lineN2line'}{$offLn}\n"; 
		}
	}
	my $stat_exon = $ms->ins_calc(\@len_exon); 
	my $stat_intron = $ms->ins_calc(\@len_intron); 
	$exonNum = $stat_exon->{'COUNT'}; 
	$avgExonLen = $stat_exon->{'interval_mean'}; 
	$avgIntronLen = $stat_intron->{'interval_mean'}; 
	( defined $avgIntronLen and $avgIntronLen ne '' ) or $avgIntronLen = -1; 
	print {$outFh} join("\t", $topID, $topS, $topE, $topLen, $exonNum, $avgExonLen, $avgIntronLen)."\n"; 
}
