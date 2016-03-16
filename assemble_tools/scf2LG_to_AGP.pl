#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"gapLen:i", # 1000
); 
$opts{'gapLen'} //= 1000; 

my $help_txt = <<HH; 

perl $0 pop1_onP1_mapI_f.csv.LGs.scf2LG pop1_onP1_mapI_f.csv.LGs.scf2LG.chr0 > pop1_onP1_mapI_f.csv.LGs.scf2LG.agp

-gapLen      [$opts{'gapLen'}]

HH
!@ARGV and &LogInforSunhh::usage($help_txt); 


my %str2char; 
%str2char = qw(
F    +
R    -
U    ?
); 

my @infor_scf2LG; 
for my $fn (@ARGV) {
	push( @infor_scf2LG, @{ &load_scf2LG_tbl($fn) } ); 
}

# [Sunhh@whale super_scaffold]$ head -3 SP_Nov_11_13_30_8_90_3_superscaffold_merged.agp
# ##agp-version   2.0
# Super_scaffold_248      1       341379  1       W       SpoScf_00639    1       341379  -
# Super_scaffold_248      341380  600696  2       N       259317  scaffold        yes     map

my %outLG; 
for my $ar1 ( @infor_scf2LG ) {
	my ($lgID, $scfID, $scfLen, $scfStr, $otherInfo) = @$ar1; 
	$outLG{'cnt'} ++; 
	my $scfStr_c = $str2char{$scfStr}; 
	defined $scfStr_c or &stopErr("[Err] Bad scfStr [$scfStr]\n"); 
	if ( defined $outLG{'LG'}{$lgID} ) {
		my @tb = @{$outLG{'LG'}{$lgID}}; 
		push( @tb, [ $outLG{'cnt'}, $lgID, $tb[-1][3]+1, $tb[-1][3]+$opts{'gapLen'},$tb[-1][4]+1,'N',$opts{'gapLen'}, 'scaffold', 'yes', 'map', 'NA' ] ); 
		push( @tb, [ $outLG{'cnt'}, $lgID, $tb[-1][3]+1, $tb[-1][3]+$scfLen,$tb[-1][4]+1,'W',$scfID, 1, $scfLen, $scfStr_c, $otherInfo ] ); 
	@{$outLG{'LG'}{$lgID}} = @tb; 
	} else {
		$outLG{'LG'}{$lgID} = [ [ $outLG{'cnt'}, $lgID, 1, $scfLen, 1, 'W', $scfID, 1, $scfLen, $scfStr_c, $otherInfo ] ]; 
	}
}
for my $lgID (sort { $outLG{'LG'}{$a}[0][0] <=> $outLG{'LG'}{$b}[0][0] } keys %{$outLG{'LG'}}) {
	my @tb = @{ $outLG{'LG'}{$lgID} }; 
	for my $tl (@tb) {
		print STDOUT join("\t", @{$tl}[ 1 .. $#$tl ])."\n"; 
	}
}

sub load_scf2LG_tbl {
	my ($fn) = @_; 
	my $fh = &openFH($fn, '<'); 
	my @back; # ( [LgID, ScfID, ScfLen, ScfDirection, ...], [], ... )
	while (<$fh>) {
		m/^\s*(#|$)/ and next; 
		chomp; s/[^\S\t]+$//; 
		my @ta = split(/\t/, "$_\n"); chomp($ta[-1]); 
		$ta[0] eq 'ChrID' and next; 
		push(@back, [ @ta[0,1,4,2,3] ]); 
	}
	close($fh); 
	return (\@back); 
}
# [Sunhh@Falcon kept]$ head pop1_onP1_mapI_f.csv.LGs.scf2LG
# ChrID   ScfID   ScfDirection    LgRange ScfLen
# Cma_Chr01       Cma_Scf00018    F       0-79.1378754030414:56   4554690
# Cma_Chr01       Cma_Scf00061    R       83.6316433435473-87.8609115785575:5     721047
# Cma_Chr01       Cma_Scf00052    U       90.8726789099144-90.8726789099144:2     947832
# Cma_Chr01       Cma_Scf00006    F       92.1737855205585-203.045084069345:79    6853530
# Cma_Chr02       Cma_Scf00041    R       0-29.8840864225909:15   1605698
# Cma_Chr02       Cma_Scf00027    R       33.3607718361044-92.6478114225589:32    3553341

