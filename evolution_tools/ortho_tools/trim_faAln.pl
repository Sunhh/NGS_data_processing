#!/usr/bin/perl
use strict; 
use warnings; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"max_gapR:f", # 1 
	"get_4d:i", 
	"as_coding!", 
	"in_fa:s", # in_aln.fas
	"out:s", # outfile 
	"wind:i", 
	"step:i", 
); 

$opts{'max_gapR'} //= 1; 
$opts{'wind'}     //= 1; 
$opts{'step'}     //= 1; 
my ($wind, $step) = ($opts{'wind'}, $opts{'step'}); 

my $help_txt = <<HH; 
################################################################################
# perl $0 -in_fa in_aln.fasta   -out trimmed_aln.fasta
#
# -help
#
# -get_4d        [translation_table_number] will set -as_coding
#
# -as_coding     [Boolean] If given, treated as codons (three bases)
#   Same to " -wind 3 -step 3 "
# -max_gapR      [$opts{'max_gapR'}] Maximum allowed gap ratio. [0-1]
################################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
defined $opts{'in_fa'} or &LogInforSunhh::usage($help_txt); 

$opts{'get_4d'} //= 0; 
my %codon_4d; 
if ( $opts{'get_4d'} > 0 ) {
	$opts{'as_coding'} = 1; 
	%codon_4d = %{ &fastaSunhh::get_4d_codon($opts{'get_4d'}) }; 
}

my $ofh_fa = \*STDOUT; 
defined $opts{'out'} and $ofh_fa = &openFH( $opts{'out'}, '>' ); 

if ($opts{'as_coding'}) {
	$wind = $step = 3; 
}

my %in_fa = %{ $fs_obj->save_seq_to_hash( 'faFile' => $opts{'in_fa'} , 'has_head'=>1) }; 
my @ids = sort { $in_fa{$a}{'Order'} <=> $in_fa{$b}{'Order'} } keys %in_fa; 
my $seqN = scalar(@ids); 
my $seqL ; 
for (@ids) { 
	$in_fa{$_}{'seq'} =~ s!\s!!g; 
	$in_fa{$_}{'len'} = length( $in_fa{$_}{'seq'} ); 
	$seqL //= $in_fa{$_}{'len'}; 
	$seqL == $in_fa{$_}{'len'} or &stopErr("[Err] Different length of sequences [$seqL VS.  $in_fa{$_}{'len'}]\n"); 
	$in_fa{$_}{'ccc'} = [ split(//, $in_fa{$_}{'seq'}) ]; 
}

my @toRM; 
for ( my $i=0; $i<$seqL; $i+=$step ) {
	my @wb; 
	for (my $k=0; $k<$seqN; $k++) {
		my $t_wb = ''; 
		for (my $j=0; $j<$wind and $j+$i<$seqL; $j++) {
			$j+$i < $seqL or last; 
			$t_wb .= $in_fa{$ids[$k]}{'ccc'}[$i+$j]; 
		}
		push(@wb, $t_wb); 
	}
	# Count for gap ratio
	my $gapN = 0; 
	for my $t_wb (@wb) {
		$t_wb =~ m![\-]! and $gapN ++; 
	}
	my $is_good_4d = 1; 
	my %frame_4d_p; 
	if ( $opts{'get_4d'} > 0 ) {
		$wind == 3 or &stopErr("[Err] wind [$wind] is not 3.\n"); 
		my %frame_4d; 
		my $cnt_total = 0; 
		for my $bbb (@wb) {
			if ( $bbb =~ m![\-]! ) {
				next; 
			}
			$cnt_total ++; 
			defined $codon_4d{$bbb} or do { $is_good_4d = 0; last; }; 
			for my $tp ( keys %{$codon_4d{$bbb}[1]} ) {
				$frame_4d{$tp} ++; 
			}
		}
		for my $tp (sort { $a <=> $b } keys %frame_4d) {
			$frame_4d{$tp} == $cnt_total and $frame_4d_p{$tp} = 1; 
		}
		scalar( keys %frame_4d_p ) > 0 or $is_good_4d = 0; 
	}
	if ( $gapN > $seqN * $opts{'max_gapR'} ) {
		for (my $j=0; $j<$wind and $j+$i<$seqL; $j++) {
			$toRM[$j+$i] //= 1; 
		}
	} elsif ( $opts{'get_4d'} > 0 ) {

		my $bbb = join('', @wb); 

		if ($is_good_4d == 1) {
			for (my $j=0; $j<$wind and $j+$i<$seqL; $j++) {
				my $tj = $j+1; 
				defined $frame_4d_p{$tj} and next; 
				$toRM[$j+$i] //= 1; 
			}
		} else {
			for (my $j=0; $j<$wind and $j+$i<$seqL; $j++) {
				$toRM[$j+$i] //= 1; 
			}
		}
	}
}

for (my $i=0; $i<$seqL; $i++) {
	defined $toRM[$i] and $toRM[$i] == 1 and next; 
	for (my $k=0; $k<$seqN; $k++) {
		$in_fa{$ids[$k]}{'got'} .= $in_fa{$ids[$k]}{'ccc'}[$i]; 
	}
}
for my $tk (@ids) {
	$in_fa{$tk}{'got'} =~ s!(.{60})!$1\n!g; 
	chomp( $in_fa{$tk}{'got'} ); 
	print {$ofh_fa} ">$tk\n$in_fa{$tk}{'got'}\n"; 
}



