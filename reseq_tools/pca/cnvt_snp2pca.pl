#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use SNP_tbl; 
my $st_obj = SNP_tbl->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"opref:s", # Default pca_in 
	"startColN:i", # Default 2 
	"insertLen:i", # 1000
); 

$opts{'opref'} //= 'pca_in'; 
$opts{'startColN'} //= 2; 
$opts{'insertLen'} //= 1000; 

my $help_txt = <<HH; 

perl $0 in.snp -opref out_prefix_for_pca_in 

-help 
-startColN      [$opts{'startColN'}]
-insertLen      [$opts{'insertLen'}]

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my @InFp ; 
-t or @InFp = ( \*STDIN ); 
if (@ARGV) {
	for (@ARGV) {
		push( @InFp, &openFH($_, '<') ); 
	}
}

my $ofh_geno = &openFH( "$opts{'opref'}.geno", '>' ); 
my $ofh_snp  = &openFH( "$opts{'opref'}.snp",  '>' ); 
my $ofh_log  = &openFH( "$opts{'opref'}.pos_log", '>' ); 
my $para_colMax; 
for my $fh (@InFp) {
	&tsmsg("[Msg] Processing [$fh]\n"); 
	my ($ref_base, $var_base); 
	my ($para_prevCID, $para_prevCEnd, $para_currCEnd); 
	$para_prevCEnd = 0; 
	LINE: 
	while (<$fh>) {
		$. % 1000 == 1 and &tsmsg("[Msg]  $. line.\n"); 
		m/^\s*(#|$)/ and next; 
		chomp; 
		my @ta = split(/\t/, $_); 
		$para_colMax //= $#ta; 
		if ($. == 1 and $ta[0] =~ m/^(chr|chrID|chromID)$/i) {
			next; 
		}
		$para_prevCID //= $ta[0]; 
		if ($para_prevCID ne $ta[0]) {
			$para_prevCID = $ta[0]; 
			$para_prevCEnd = $para_currCEnd + $opts{'insertLen'}; 
		}
		$para_currCEnd = $ta[1]; 

		print {$ofh_snp} "$ta[0]_MRK_$ta[1]\t1\t0.0\t" . ($ta[1]+$para_prevCEnd) ; 
		print {$ofh_log} join("\t", "$ta[0]_MRK_$ta[1]", $ta[0], $ta[1], 1, $ta[1]+$para_prevCEnd)."\n"; 

		($ref_base, $var_base) = ('', ''); 
		my @tmp_al; 
		for (my $i=$opts{'startColN'}; $i<=$para_colMax; $i++) {
			$ta[$i] //= 'N'; 
			length( $ta[$i] ) > 0 or $ta[$i] = 'N'; 
			$ta[$i] eq '-' and $ta[$i] = 'N'; 
			$ta[$i] eq 'N' and do { next; }; 
			my $tb = $st_obj->SingleChar( $ta[$i], 'maxAlleleN'=>2, 'onlyATGC'=>0 ) ; 
			@{$tmp_al[$i]} = &SNP_tbl::dna_d2b($tb); 
			$tmp_al[$i][1] //= $tmp_al[$i][0]; 
			unless ($ref_base ne '' and $var_base ne '') {
				for my $tc (@{$tmp_al[$i]}) {
					$ref_base eq '' and $ref_base = $tc; 
					$ref_base ne $tc and do { $var_base = $tc; last; }; 
				}
			}
		}
		print {$ofh_snp} "\t$ref_base\t$var_base\n"; 
		for (my $i=$opts{'startColN'}; $i<=$para_colMax; $i++) {
			$ta[$i] eq 'N' and do { print {$ofh_geno} '9'; next; }; 
			my $cnt = 0; 
			for my $tb ( @{$tmp_al[$i]}[0,1] ) {
				$tb eq $ref_base and $cnt++; 
			}
			print {$ofh_geno} $cnt;
		}
		print {$ofh_geno} "\n"; 
	}
}
close($ofh_geno); 
close($ofh_snp); 


