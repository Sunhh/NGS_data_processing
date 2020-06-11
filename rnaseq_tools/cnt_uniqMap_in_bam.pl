#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use SeqAlnSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"inBam:s@", 
	"exe_samtools:s", 
); 

my %gg; 
$gg{'help_txt'} = <<"H1"; 
################################################################################
# perl $0 -inBam in.bam > in.bam.uniq_multi_hits_number
#
# -help
#
# -exe_samtools       ['samtools']
#
#######
# Result file (in.bam.uniq_multi_hits_number) example: 
# InBam                    RS_all  RS_unique    RS_multiple   RS_unmapped    RP_all       RP_unique       RP_multiple    RP_unmapped
# F582R1_bySM_fix.bam      0       0            0             0              28913649     28740150        123542         49957
# F582R2HG_bySM_fix.bam    0       0            0             0              28241020     27822134        334698         84188
# SQ026_bySM_fix.bam       0       0            0             0              142786057    140123579       833960         1828518
################################################################################
H1
$gg{'exe_samtools'} = 'samtools'; 
defined $opts{'exe_samtools'} and $gg{'exe_samtools'} = $opts{'exe_samtools'}; 
$opts{'help'} and &LogInforSunhh::usage($gg{'help_txt'}); 
defined $opts{'inBam'} or &LogInforSunhh::usage($gg{'help_txt'}); 


my %flag_R1 = %{ &SeqAlnSunhh::mk_flag( 'keep' => '0=1,6=1,7=0' ) }; 
my %flag_R2 = %{ &SeqAlnSunhh::mk_flag( 'keep' => '0=1,6=0,7=1' ) }; 
my %flag_RS = %{ &SeqAlnSunhh::mk_flag( 'keep' => '0=0' ) }; 
my %flag_UN = %{ &SeqAlnSunhh::mk_flag( 'keep' => '2=1' ) };
my %flag2R; 
for (keys %flag_R1) { $flag2R{$_} = 'R1'; }
for (keys %flag_R2) { $flag2R{$_} = 'R2'; }
for (keys %flag_RS) { $flag2R{$_} = 'RS'; }

my %rdInfo; 
for my $ibam (@{$opts{'inBam'}}) {
	my $fh = &SeqAlnSunhh::openSam( $ibam, undef(), { 'wiH'=>0, 'exe_samtools'=>$gg{'exe_samtools'} } ); 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		defined $flag2R{$ta[1]} or &stopErr("[Err] Failed to understand flag [$ta[1]] in line: $_\n"); 
		if ( defined $flag_UN{$ta[1]} ) {
			$rdInfo{$ibam}{$ta[0]}{$flag2R{$ta[1]}} = 0; 
		} else {
			$rdInfo{$ibam}{$ta[0]}{$flag2R{$ta[1]}} ++; 
		}
	}
	close($fh); 
}
my %cnt; 
for my $ibam (keys %rdInfo) {
	for my $rdID (keys %{$rdInfo{$ibam}}) {
		if ((defined $rdInfo{$ibam}{$rdID}{'R1'} or defined $rdInfo{$ibam}{$rdID}{'R2'}) and defined $rdInfo{$ibam}{$rdID}{'RS'}) {
			&stopErr("[Err] rdID [$rdID] has R1/R2 and RS\n"); 
		} elsif ( defined $rdInfo{$ibam}{$rdID}{'RS'} ) {
			if ( $rdInfo{$ibam}{$rdID}{'RS'} == 1 ) {
				$cnt{'RS_unique'} ++; 
			} elsif ($rdInfo{$ibam}{$rdID}{'RS'} == 0) {
				$cnt{'RS_unmapped'} ++; 
			} else {
				$cnt{'RS_multiple'} ++; 
			}
		} else {
			if (defined $rdInfo{$ibam}{$rdID}{'R1'} and defined $rdInfo{$ibam}{$rdID}{'R2'}) {
				if ( $rdInfo{$ibam}{$rdID}{'R1'} > 1 or $rdInfo{$ibam}{$rdID}{'R2'} > 1 ) {
					$cnt{'RP_multiple'} ++; 
				} elsif ( $rdInfo{$ibam}{$rdID}{'R1'} == 0 and $rdInfo{$ibam}{$rdID}{'R2'} == 0 ) {
					$cnt{'RP_unmapped'} ++; 
				} else {
					$cnt{'RP_unique'} ++; 
				}
			} elsif (defined $rdInfo{$ibam}{$rdID}{'R1'}) {
				if ( $rdInfo{$ibam}{$rdID}{'R1'} == 1 ) {
					$cnt{'RP_unique'} ++; 
				} elsif ( $rdInfo{$ibam}{$rdID}{'R1'} == 0 ) {
					$cnt{'RP_unmapped'} ++; 
				} else {
					$cnt{'RP_multiple'} ++; 
				}
			} elsif (defined $rdInfo{$ibam}{$rdID}{'R2'}) {
				if ( $rdInfo{$ibam}{$rdID}{'R2'} == 1 ) {
					$cnt{'RP_unique'} ++; 
				} elsif ( $rdInfo{$ibam}{$rdID}{'R2'} == 0 ) {
					$cnt{'RP_unmapped'} ++;
				} else {
					$cnt{'RP_multiple'} ++; 
				}
			} else {
				&stopErr("[Err] Why here?\n"); 
			}
		}
	}
}

for (qw/RS_unique RS_multiple RS_unmapped RP_unique RP_multiple RP_unmapped/) {
	$cnt{$_} //= 0; 
}

$cnt{'RS_all'} = $cnt{'RS_unique'} + $cnt{'RS_multiple'} + $cnt{'RS_unmapped'}; 
$cnt{'RP_all'} = $cnt{'RP_unique'} + $cnt{'RP_multiple'} + $cnt{'RP_unmapped'}; 

print join("\t", qw/InBam RS_all RS_unique RS_multiple RS_unmapped RP_all RP_unique RP_multiple RP_unmapped/)."\n"; 
print join("\t", join(";", @{$opts{'inBam'}}), @cnt{qw/RS_all RS_unique RS_multiple RS_unmapped RP_all RP_unique RP_multiple RP_unmapped/})."\n"; 

