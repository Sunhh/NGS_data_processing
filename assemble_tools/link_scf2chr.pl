#!/usr/bin/perl 
use strict; 
use warnings; 
!@ARGV and die "perl $0 inFa inLis\n"; 
# CG_ChrID        Final_ID        Final_Direction Ordered Annotation
# CG_Chr01        CGscf_0024      R       Y
# CG_Chr01        CGscf_0019      R       Y

my $faF = shift; 
my $lisF = shift; 

my $preTag='Mos'; 

my %seq; 
my %len; 
my $id; 
open FF,'<',"$faF" or die; 
while (<FF>) {
	chomp; 
	if (m/^>(\S+)/) {
		$id = $1; 
		defined $seq{$id} and die "repeat=$id\n"; 
	}else{
		$seq{$id} .= "$_"; 
	}
}
close FF; 

for (keys %seq) {
	$seq{$_} =~ s/\s//g; 
	$len{$_} = length( $seq{$_} ); 
}


my $split_N = 'N' x 1e3; 

open LF,'<',"$lisF" or die; 
my %used; 
my %chr_seq; 
my %chr_len; 
while (<LF>) {
	chomp; /^\s*$/ and next; 
	my @ta = split(/\t/, $_); 
	($ta[0] eq 'CG_ChrID' or $ta[0] eq 'ChrID') and next; 
	defined $seq{$ta[1]} or die "?$ta[1]?\n"; 
	my $add_seq = $seq{$ta[1]}; 
	if ( $ta[2] eq 'F' or $ta[2] eq 'U' or $ta[2] eq 'N' ) {
	} elsif ( $ta[2] eq 'R' ) {
		$add_seq = reverse($add_seq); 
		$add_seq =~ tr/ATGCatgc/TACGtacg/; 
	} else {
		die "t1: $_\n"; 
	}
	$used{$ta[1]} = 1; 

	if (defined $chr_seq{$ta[0]}) {
		# print STDERR join("\t", $ta[0], $ta[1], length($chr_seq{$ta[0]})+1, length($chr_seq{$ta[0]})+length($add_seq)+length($split_N), $ta[2], $ta[3], length($add_seq) )."\n"; 
		print STDERR join("\t", $ta[0], $ta[1], length($chr_seq{$ta[0]})+1+length($split_N), length($chr_seq{$ta[0]})+length($add_seq)+length($split_N), $ta[2], $ta[3], length($add_seq) )."\n"; 
		$chr_seq{$ta[0]} .= "$split_N$add_seq"; 
	} else {
		print STDERR join("\t", $ta[0], $ta[1],                          1, length($add_seq), $ta[2], $ta[3], length($add_seq))."\n"; 
		$chr_seq{$ta[0]} = $add_seq; 
	}
}
close LF; 

my $rest; 
for my $td ( sort { $len{$b} <=> $len{$a} || $a cmp $b } keys %seq ) {
	if (!defined $used{$td}) {
		if (defined $rest and $rest ne '') {
			print STDERR join("\t", "${preTag}_Chr00", $td, length($rest)+1, length($rest)+length($split_N)+length($seq{$td}), qw/U N/, length($seq{$td}) )."\n"; 
			$rest .= "$split_N$seq{$td}"; 
		}else{
			print STDERR join("\t", "${preTag}_Chr00", $td, 1, length($seq{$td}), qw/U N/, length($seq{$td}))."\n"; 
			$rest = $seq{$td}; 
		}
	}
}

for my $td (sort keys %chr_seq) {
	my $oseq = $chr_seq{$td}; 
	$oseq =~ s/(.{100})/$1\n/g; chomp($oseq);
	print STDOUT ">$td\n$oseq\n"; 
}
$rest =~ s/(.{100})/$1\n/g; chomp($rest); 
print STDOUT ">${preTag}_Chr00\n$rest\n"; 



