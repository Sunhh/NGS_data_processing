#!/usr/bin/perl 
use strict; 

# jellyfish-2.0 query P1_Cor1st_1k_m81.jf GTCGATTAAATACGCAAGTTTGCCACTAGCTGTAGCAGCACTGTCTACGAGTCACGAGATGATGTACAAATGACCGCATAGCGCAGTGGATTAGCGCCTCTGACTTCGGATCAGAAGGTTGTGGGTTTGACTCCCACTGTGGTCATCTGT
my $db_jf = 'P1_Cor1st_1k_m61.jf'; 
$db_jf = 'P1_Cor1st_3h5h1k_m61.jf'; 
my $db_km = 61; 
my $exe_jf2 = '/home/Sunhh/src/Assemble/MaSuRCA/MaSuRCA-2.2.2/bin/jellyfish-2.0'; 

my $step = 10; 
$step = 1; 

my (@IDs, %seqs); 
my $tk='unknown'; 
if (!(-t)) {
	while (<>) {
		if (m/^>(\S+)/) {
			$tk = $1; 
			$tk eq 'unknown' and die "Should no be unknown\n"; 
			defined $seqs{$tk} and die "repeat $tk\n"; 
		}else{
			defined $seqs{$tk} or push(@IDs, $tk); 
			$seqs{$tk} .= $_; 
		}
	}
}else{
	for my $fn (@ARGV) {
		if (-f $fn) {
			open F,'<',"$fn" or die; 
			while (<F>) {
				if (m/^>(\S+)/) {
					$tk = $1; 
					$tk eq 'unknown' and die "Should no be unknown\n"; 
					defined $seqs{$tk} and die "repeat $tk\n"; 
				}else{
					defined $seqs{$tk} or push(@IDs, $tk); 
					$seqs{$tk} .= $_; 
				}
			}
			close F; 
		}else{
			defined $seqs{$tk} or push(@IDs, $tk); 
			$seqs{$tk} .= $fn; 
		}
	}
}

# defined $seqs{'unknown'} and unshift(@IDs, 'unknown'); 
for (@IDs) {
	$seqs{$_} =~ s/\s//g; 
}
for my $tk (@IDs) {
print STDOUT ">$tk\n"; 
	my $len = length($seqs{$tk}); 
	if ($len > 0) {
		for (my $i=0; $i+$db_km<=$len; $i+=$step) {
			my $kseq = substr($seqs{$tk}, $i, $db_km); 
			my $qry_c = `$exe_jf2 query $db_jf $kseq`; 
			$qry_c =~ m/^(\S+)\s+(\d+)$/ or die "fasd\n"; 
			my ($res_seq, $res_c) = ($1, $2); 
			if ( $res_seq eq $kseq ) {
			} else {
				my $rc_kseq = &rc($kseq); 
				$rc_kseq eq $res_seq or $res_c = -1; 
			}
			print STDOUT join("\t", $i+1, $kseq, $res_c)."\n"; 
		}
	}
}

sub rc {
	$_ = shift; 
	$_ =~ tr/ATGCatgc/TACGtacg/; 
	return reverse($_); 
}


sub tsmsg {
	my $tt = scalar(localtime()); 
	print STDERR join('', "[$tt]", @_); 
}
