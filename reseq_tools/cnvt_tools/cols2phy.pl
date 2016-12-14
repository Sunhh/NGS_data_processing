#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
use SNP_tbl; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"colN_start:i", # 2 
	"input_is_tab!", 
); 
$opts{'colN_start'} //= 2; 

my $help_txt = <<HH; 
######################################################################
# perl $0 in.cols.snp > in.cols.snp.phy
#
#   -input_is_tab    [Bool] 
#   -colN_start      [$opts{'colN_start'}]
#
######################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

my $colS = $opts{'colN_start'}; 
my %d2b_list;
my %b2d_list; 
{
my @aa = (
[qw/W A T/],
[qw/S C G/],
[qw/M A C/],
[qw/K G T/],
[qw/R A G/],
[qw/Y C T/],
[qw/A A A/], 
[qw/T T T/], 
[qw/G G G/], 
[qw/C C C/], 
);
for my $tr (@aa) {
	my @bb = @$tr;
	$d2b_list{$bb[0]} = [$bb[1], $bb[2]];
	$b2d_list{"$bb[1]$bb[2]"} = $bb[0]; 
	$b2d_list{"$bb[2]$bb[1]"} = $bb[0]; 
}
}

my %used; 

my @header; 
my @seqs; 
SNP_LINE: 
while (<>) {
	s/[^\S\t]+$//; 
	my @ta = split(/\t/, $_); 
	if ( $. == 1 ) {
		@header=@ta; 
		next; 
	}
	for (my $i=$colS; $i<@ta; $i++) {
		$ta[$i] = uc($ta[$i]); 
		if ( $opts{'input_is_tab'} ) {
			if ( $ta[$i] =~ m!^([ATGC])/([ATGC])$! ) {
				$ta[$i] = $b2d_list{"$1$2"}; 
			} elsif ( $ta[$i] eq './.' or $ta[$i] eq 'N/N' ) {
				$ta[$i] = 'N'; 
			} elsif ( $ta[$i] =~ m!^([ATGCN*]+)/([ATGCN*]+)$! ) {
				&tsmsg("[Wrn] Skip site [$ta[0] $ta[1]] with indel [$ta[$i]]\n"); 
				next SNP_LINE; 
			} else {
				&stopErr("[Err] Bad genotype [$ta[$i]]\n"); 
			}
		} elsif ( $ta[$i] =~ m!^[ATGC]$! ) {
			; 
		} elsif ( $ta[$i] =~ m!^([ATGC])([ATGC])$! ) {
			$ta[$i] = $b2d_list{$ta[$i]}; 
		} elsif ( $ta[$i] =~ m![+*]! ) {
			&tsmsg("[Wrn] Skip site [$ta[0] $ta[1]] with indel [$ta[$i]]\n"); 
			next SNP_LINE; 
		} else {
			unless ( defined $used{'bad_geno'}{$ta[$i]} ) {
				&tsmsg("[Wrn] Skip SNP site with bad genotype [$ta[$i]]\n"); 
			}
			&tsmsg("[Wrn] Skip site [$ta[0] $ta[1]] with bad genotype [$ta[$i]]\n"); 
			next SNP_LINE; 
		}
	}
	for (my $i=$colS; $i<@ta; $i++) {
		$seqs[$i] .= $ta[$i]; 
	}
}

for (@seqs) {
	defined $_ or next; 
	$_ =~ s!\s!!g; 
}
my $indv_N = $#header-$colS+1; 
my $base_N = length($seqs[$colS]); 

print STDOUT <<OOO; 
   $indv_N   $base_N
OOO

my %used_id; 
my %has_out_i; 
for (my $p=0; $p<$base_N; $p+=50) {
	for (my $i=$colS; $i<@header; $i++) {
		my $indv_ID = ''; 
		unless ( defined $has_out_i{$i} ) {
			$has_out_i{$i} = 1; 
			$indv_ID = $header[$i]; 
			$indv_ID =~ s!\s!_!g; 
			$indv_ID =~ tr!)(][:;,!_______!; 
			length($indv_ID) > 9 and $indv_ID = substr($indv_ID, 0, 9); 
			$indv_ID = sprintf("%-10s", $indv_ID); 
		        defined $used_id{$indv_ID} and die "Repeat ID [$indv_ID]\n"; 
			$used_id{$indv_ID} = 1; 
		}
		print STDOUT "${indv_ID}"; 
		my @o_seq; 
		for (my $s=$p; $s<$p+50 and $s < $base_N; $s+=10) {
			my $e = $s+10; 
			$e <= $base_N or $e = $base_N; 
			push(@o_seq, substr($seqs[$i], $s, $e-$s)); 
		}
		print STDOUT join(" ", @o_seq)."\n"; 
	}
	print STDOUT "\n"; 
}


