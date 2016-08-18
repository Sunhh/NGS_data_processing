#!/usr/bin/perl 
use strict; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"jf_db:s", # Required. 
	"jf_klen:i", # Default 17
	"jf_exe:s", # Default jellyfish
	"step:i", # Default 1 
); 

&_prepare_para(); 

sub _prepare_para {
	$opts{'jf_exe'}  //= 'jellyfish'; 
	$opts{'jf_klen'} //= 17; 
	$opts{'step'}    //= 1; 
}

my $help_txt = <<HH; 

perl $0   -jf_db  in.jf_database 

-jf_exe       [$opts{'jf_exe'}]
-jf_klen      [$opts{'jf_klen'}]
-step         [$opts{'step'}]

HH

defined $opts{'jf_db'} or &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt);


my $db_jf = $opts{'jf_db'}; 
my $db_km = $opts{'jf_klen'}; 
my $exe_jf2 = $opts{'jf_exe'}; 
my $step = $opts{'step'}; 


# jellyfish-2.0 query P1_Cor1st_1k_m81.jf GTCGATTAAATACGCAAGTTTGCCACTAGCTGTAGCAGCACTGTCTACGAGTCACGAGATGATGTACAAATGACCGCATAGCGCAGTGGATTAGCGCCTCTGACTTCGGATCAGAAGGTTGTGGGTTTGACTCCCACTGTGGTCATCTGT

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
			if ( uc($res_seq) eq uc($kseq) ) {
			} else {
				my $rc_kseq = &rc($kseq); 
				uc($rc_kseq) eq uc($res_seq) or $res_c = -1; 
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
