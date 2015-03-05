#!/usr/bin/perl -w 
# Edit at 20130117: '*' stands for a deletion, and the reads supporting deletion are calculated. 
use strict; 
# WM97_Chr01      443     G       7       .,,.,,. AJCJJJI
# WM97_Chr01      444     T       7       .,,.,,. ;ICFJJG
# WM97_Chr01      445     A       8       .$,,$.,,.^].    @JCHJJIB
# WM97_Chr01      446     T       6       ,.,,..  JIJJJC

-t and !@ARGV and die "perl $0 input.pileup\n"; 

my $refBaseC  =   2; 
my $readBaseC =   4; 
my $covC      =   3; 
my $qualC     =   5;
my $qualCut   = -99; 
my $covCut    =   1; 

my @genotype  = qw/A T G C N/; 

# print STDOUT join("\t", qw/RefChr ChrPos RefBase RdCov/, @genotype, 'TotalUsed', 'TotalDrop', 'MinorInMaxAndMinor', 'Max2RatByAlphabet')."\n"; 
print STDOUT join("\t", qw/RefChr ChrPos RefBase RdCov/, @genotype, 'Deletion', 'MinorInMaxAndMinor', 'Max2RatByAlphabet')."\n"; 
while (<>) {
	$. % 1000000 == 1 and warn "[Stat]$. lines.\n"; 
	chomp; m/^\s*$/ and next; m/^\#/ and next; 
	my @ta = split(/\t/, $_); 
	if ($ta[$covC] < $covCut) {
		print STDOUT join("\t", @ta[0,1,$refBaseC,$covC], qw/0 0 0 0 0 0 - -/, $ta[$covC])."\n"; 
	} else {
		my $ref_base   = uc($ta[$refBaseC]); 
		my $read_bases = $ta[$readBaseC]; 
		my $read_quals = $ta[$qualC]; 
		if ($read_bases =~ m/[\$\^\+-]/) {
			$read_bases =~ s/\^.//g; #removing the start of the read segement mark
			$read_bases =~ s/\$//g; #removing end of the read segment mark
			while ($read_bases =~ m/[\+-](\d+)/g) {
				my $indel_len = $1;
				$read_bases =~ s/[\+-]$indel_len.{$indel_len}//; # remove indel info from read base field
			}
		}
		if ( length($read_bases) != length($read_quals) ) {
			warn "[Error]Skip line: $_\n"; 
		}
		# after removing read block and indel data the length of read_base 
		# field should identical to the length of read_quals field
		my @bases = split(//, $read_bases);
		my @quals; 
		if ($qualCut > 0) {
			@quals  = split(//, $read_quals); 
			my @new_bases; 
			my @new_quals; 
			for (my $i=0; $i<@bases; $i++) {
				if ( ord($quals[$i])-33 >= $qualCut and $bases[$i] ne '*' ) {
					# '*' here means a read supporting deletion. 
					if ($bases[$i] =~ m/[ATGCN]/i) {
						push(@new_bases, $bases[$i]); 
					} elsif ($bases[$i] =~ m/[.,]/) {
						push(@new_bases, $ref_base);
					}
					push(@new_quals, $quals[$i]); 
				}
			}
			@quals = @new_quals; 
			@bases = @new_bases; 
		}
		my %count; 
		for (my $i=0; $i<@bases; $i++) {
			if ($bases[$i] =~ m/[.,]/) {
				$count{ $ref_base } ++; 
			}else{
				$count{ uc($bases[$i]) }++; 
			}
		}
		my $ttl_used = 0; 
		for my $tbase (@genotype) {
			defined $count{$tbase} or $count{$tbase} = 0; 
			$ttl_used += $count{$tbase}; 
		}
		my $deletion_num = $ta[$covC]-$ttl_used; 
		$count{'deletion'} = $deletion_num; 
		my @tbs1 = sort { $count{$b} <=> $count{$a} || $a cmp $b } grep { $_ ne 'N' } (@genotype, 'deletion'); 
		my @tbs2 = sort { $a cmp $b } @tbs1[0,1]; 
		my ($rat1, $rat2) = ('-', '-'); 
		if ($count{$tbs1[0]} > 0) {
			my $bi_sum = $count{$tbs1[0]}+$count{$tbs1[1]}; 
			$rat1 = $count{$tbs1[1]}/$bi_sum; 
			$rat2 = $count{$tbs2[0]}/$bi_sum; 
		}
#		print STDOUT join("\t", @ta[0,1,$refBaseC,$covC], @count{@genotype}, $ttl_used, $ta[$covC]-$ttl_used, $rat1, $rat2)."\n"; 
		print STDOUT join("\t", @ta[0,1,$refBaseC,$covC], @count{@genotype}, $deletion_num, $rat1, $rat2)."\n"; 
	}
}



