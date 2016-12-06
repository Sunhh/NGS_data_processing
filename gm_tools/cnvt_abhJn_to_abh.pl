#!/usr/bin/perl
### Format of .abh.jn : 
# Offs_ID ChrID   StartP  EndP    Offs_abh        Offs_siteN2     Pop_siteN1      Pop_siteN2
# RIL1    chr00   1912    10453519        a       910     1213    910
# RIL1    chr00   10623549        12245749        b       76      106     76
# RIL1    chr00   13142934        15068477        a       66      96      66
# RIL1    chr00   15191407        16135441        b       63      71      63
# RIL1    chr00   17209079        18823573        b       254     373     254
# RIL1    chr00   18935371        20542313        a       81      101     81
### Format of .abh : 
# chr     pos     base    Dulce   PI414723        RIL1    RIL2    RIL3    RIL4    RIL5    RIL6    RIL7    RIL9    RIL10   RIL11   RIL12   RIL13   RIL14   RIL15   RIL16
# chr00   1912    T       T/T     C/C     a       b       b       a       a       a       a       a       a       b       b       h       h       h       a       a
# chr00   2002    T       T/T     A/A     a       b       b       a       a       a       a       a       a       b       b       h       h       h       a       a
# chr00   2036    A       A/A     G/G     a       b       b       a       a       a       a       a       a       b       b       h       a       h       a       a
# chr00   2074    C       C/C     G/G     a       b       b       a       a       a       a       a       a       b       b       h       a       b       a       a
# chr00   4267    T       T/T     A/A     a       b       b       a       a       a       a       a       a       b       b       h       a       b       a       a
# chr00   4316    G       G/G     A/A     a       b       b       a       a       a       a       a       a       b       b       h       a       b       a       a
# chr00   4591    C       C/C     A/A     a       b       b       a       a       a       a       a       a       b       b       h       a       b       a       a
# chr00   5238    T       T/T     C/C     a       b       b       a       a       a       a       a       a       b       b       h       a       b       a       a
# chr00   6342    C       C/C     T/T     a       b       b       a       a       a       a       a       a       b       b       h       a       b       a       a
# chr00   6559    A       A/A     C/C     a       b       b       a       a       a       a       a       a       b       b       h       a       b       a       a
# chr00   11011   C       C/C     T/T     a       b       b       a       a       a       a       a       a       b       b       h       a       b       a       a
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
!@ARGV and die "perl $0 in.abh.jn [in.abh] > in.abh.jn.abh\n"; 

my $jnF = shift; 
my $src_abhF = shift; 
my %neighbor; 
if ( defined $src_abhF ) {
	%neighbor = &_load_abh_pos( $src_abhF ); 
}
sub _load_abh_pos {
	my $fn = shift; 
	my $fh = &openFH($fn, '<'); 
	my @prev; 
	my %back; 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		if ( @prev > 0 and $prev[0] eq $ta[0]) {
			my $tk = "$prev[1]\t$ta[1]"; 
			$back{$prev[0]}{$tk} = 1; 
		}
		@prev = @ta[0,1]; 
	}
	close ($fh); 
	return(%back); 
}# _load_abh_pos () 


my $jnH = &openFH($jnF, '<'); 
my %info; 
my (%loci, %lociH); 
my @indv; 
while (<$jnH>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[0] eq 'Offs_ID' and next; 
	defined $info{$ta[0]} or push(@indv, $ta[0]); 
	push(@{$info{$ta[0]}{$ta[1]}}, [@ta[2..$#ta]]); 
	$lociH{$ta[1]}{$ta[2]} = 1; 
	$lociH{$ta[1]}{$ta[3]} = 1; 
}
close($jnH); 
#### Make intervals 
for my $cid (sort keys %lociH) {
	&tsmsg("[Msg] Sorting [$cid]\n"); 
	@{$loci{$cid}} = sort { $a <=> $b } keys %{$lociH{$cid}}; 
}
#### Fill loci
print STDOUT join("\t", qw/chr pos len Spos Epos/, @indv)."\n"; 
for my $cid (sort keys %loci) {
	my @pp = @{$loci{$cid}}; 
	&tsmsg("[Msg] Filling [$cid] for [$#pp] windows\n"); 
	for (my $i=0; $i<@pp; $i++) {
		my @se_pairs; 
		push(@se_pairs, [$pp[$i]-0.5, $pp[$i]+0.5]); 
		if ( $i < $#pp ) {
			my $tk = "$pp[$i]\t$pp[$i+1]"; 
			unless ( defined $neighbor{$cid}{$tk} ) {
				push(@se_pairs, [$pp[$i]+0.5, $pp[$i+1]-0.5]); 
			}
			
		}
		for my $t1 ( @se_pairs ) {
			my ($s, $e) = @$t1; 
			# $e-$s > 100 or next; 
			my @geno; 
			my $is_only_u = 1; 
			INDV: 
			for my $iid ( @indv ) {
				for my $tr ( @{$info{$iid}{$cid}} ) {
					$tr->[0] > $e and last; 
					$tr->[1] < $s and next; 
					( $tr->[0] <= $s and $tr->[1] >= $e ) or $tr->[0]-0.5 == $s or $tr->[1]+0.5 == $e or &stopErr("[Err] Bad region [$s-$e] in [$tr->[0] $tr->[1]]\n"); 
					push(@geno, $tr->[2]); 
					$is_only_u = 0; 
					next INDV; 
				}
				# &tsmsg( "Bad: $iid [$cid $s $e]\n" ); 
				push(@geno, 'u'); 
			}
			$is_only_u or print STDOUT join("\t", $cid, $s+0.5, $e-$s, $s+0.5, $e-0.5, @geno)."\n"; 
		}
	}
}





