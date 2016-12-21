#!/usr/bin/perl
# 2016-12-15 add .vcf.tab format support, and only accept bi-allele line. 
# InDel is not supported yet. 
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
	"input_is_tab!", 
	"use_abnormal_line!", 
); 

$opts{'opref'} //= 'pca_in'; 
$opts{'startColN'} //= 2; 
$opts{'insertLen'} //= 1000; 

my $help_txt = <<HH; 

perl $0 in.snp -opref out_prefix_for_pca_in 

-help 
-startColN          [$opts{'startColN'}]
-insertLen          [$opts{'insertLen'}]

-input_is_tab       [Boolean] Required if input file is .vcf.tab format; 
-use_abnormal_line  [Boolean] If not given, site lines with InDel or abnormal genotypes are ignored. 

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my %d2b_list = &SNP_tbl::get_diploid_d2b(); 
my %used; 


my @InFp ; 
-t or @InFp = ( \*STDIN ); 
if (@ARGV) {
	for (@ARGV) {
		push( @InFp, &openFH($_, '<') ); 
	}
}
!@InFp and &LogInforSunhh::usage($help_txt); 

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
			next LINE; 
		}
		$para_prevCID //= $ta[0]; 
		if ($para_prevCID ne $ta[0]) {
			$para_prevCID = $ta[0]; 
			$para_prevCEnd += ($para_currCEnd + $opts{'insertLen'}); 
		}
		$para_currCEnd = $ta[1]; 

		($ref_base, $var_base) = ('', ''); 
		my @tmp_al; 
		for (my $i=$opts{'startColN'}; $i<=$para_colMax; $i++) {
			$ta[$i] = uc($ta[$i]); 
			if ( $opts{'input_is_tab'} ) {
				$ta[$i] //= './.'; 
				length($ta[$i]) > 0 or $ta[$i] = './.'; 
				$ta[$i] =~ m!^([^\s/]+)/([^\s/]+)$! or &stopErr("[Err] Bad genotype 01 [$ta[$i]]\n"); 
				my ($a1, $a2) = ($1, $2); 
				$a1 eq '.' and next; 
				if ( $a1 =~ m!^[ATGC]$! and $a2 =~ m!^[ATGC]$! ) {
					@{$tmp_al[$i]} = ($a1, $a2); 
					unless ($ref_base ne '' and $var_base ne '') {
						for my $tc (@{$tmp_al[$i]}) {
							$ref_base eq '' and $ref_base = $tc; 
							$ref_base ne $tc and do { $var_base = $tc; last; }; 
						}
					}
				} elsif ( $a1 =~ m!^[ATGCN*\-]+$! and $a2 =~ m!^[ATGCN*\-]+$! ) {
					if ( $opts{'use_abnormal_line'} ) {
						unless ( defined $used{'bad_geno'}{$ta[$i]} ) {
							&tsmsg("[Wrn]   Mask bad genotype [$ta[$i]]\n"); 
							$used{'bad_geno'}{$ta[$i]} = 1; 
						}
						$ta[$i] = './.'; 
					} else {
						next LINE; 
					}
				} else {
					&stopErr("[Err]   Unknown genotype 02 [$ta[$i]]\n"); 
				}
			} else {
				$ta[$i] //= 'N'; 
				length( $ta[$i] ) > 0 or $ta[$i] = 'N'; 
				$ta[$i] eq '-' and $ta[$i] = 'N'; 
				$ta[$i] eq 'N' and do { next; }; 
				if ( $ta[$i] =~ m!^[ATGC]$! ) {
					@{$tmp_al[$i]} = ($ta[$i], $ta[$i]); 
				} elsif ( $ta[$i] =~ m!^([ATGC])([ATGC])$! ) {
					@{$tmp_al[$i]} = ($1, $2); 
				} elsif ( defined $d2b_list{$ta[$i]} ) {
					@{$tmp_al[$i]} = @{ $d2b_list{$ta[$i]} }; 
				} elsif ( $opts{'use_abnormal_line'} ) {
					unless ( defined $used{'bad_geno'}{$ta[$i]} ) {
						&tsmsg("[Wrn]   Mask bad genotype 03 [$ta[$i]]\n"); 
						$used{'bad_geno'}{$ta[$i]} = 1; 
					}
					$ta[$i] = 'N'; 
				} else {
					next LINE; 
				}
				@{$tmp_al[$i]} > 0 or next; 
				unless ($ref_base ne '' and $var_base ne '') {
					for my $tc (@{$tmp_al[$i]}) {
						$ref_base eq '' and $ref_base = $tc; 
						$ref_base ne $tc and do { $var_base = $tc; last; }; 
					}
				}
			}
			if ( !($opts{'use_abnormal_line'}) and $ref_base ne '' and $var_base ne '' and @{$tmp_al[$i]} > 0 ) {
				# There are more than two alleles. 
				for my $tc (@{$tmp_al[$i]}) {
					$tc eq $ref_base and next; 
					$tc eq $var_base and next; 
					next LINE; 
				}
			}
		}
		( defined $var_base and $var_base ne '' ) or next LINE; 
		print {$ofh_snp} "$ta[0]_MRK_$ta[1]\t1\t0.0\t" . ($ta[1]+$para_prevCEnd) ; 
		print {$ofh_log} join("\t", "$ta[0]_MRK_$ta[1]", $ta[0], $ta[1], 1, $ta[1]+$para_prevCEnd)."\n"; 

		print {$ofh_snp} "\t$ref_base\t$var_base\n"; 
		for (my $i=$opts{'startColN'}; $i<=$para_colMax; $i++) {
			if ( $opts{'input_is_tab'} ) {
				$ta[$i] eq './.' and do { print {$ofh_geno} '9'; next; }; 
			} else {
				$ta[$i] eq 'N' and do { print {$ofh_geno} '9'; next; }; 
			}
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


