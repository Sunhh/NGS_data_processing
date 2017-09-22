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
	"chr_keyLen:s", # file from "deal_fasta.pl -attr key:len\n"
	"cpuN:i", 
); 

$opts{'opref'} //= 'pca_in'; 
$opts{'startColN'} //= 2; 
$opts{'insertLen'} //= 1000; 
$opts{'cpuN'}      //= 20; 

my $help_txt = <<HH; 

perl $0 in.snp -opref out_prefix_for_pca_in 

-help 
-startColN          [$opts{'startColN'}]
-insertLen          [$opts{'insertLen'}]

-input_is_tab       [Boolean] Required if input file is .vcf.tab format; 
-use_abnormal_line  [Boolean] If not given, site lines with InDel or abnormal genotypes are ignored. 

-chr_keyLen         [filename] (Optional.) Format: chr_key \\t chr_length \\n

-cpuN               [$opts{'cpuN'}] Not used yet. 

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my %d2b_list = &SNP_tbl::get_diploid_d2b(); 
my %used; 


my @InFp ; 
-t or @InFp = ( \*STDIN ); 
if (@ARGV) {
	@InFp = (); 
	for (@ARGV) {
		push( @InFp, &openFH($_, '<') ); 
	}
}
!@InFp and &LogInforSunhh::usage($help_txt); 


my $pm = &LogInforSunhh::get_pm( $opts{'cpuN'} ); 

my %chr_start; 
if (defined $opts{'chr_keyLen'}) {
	my $prevE = 0; 
	for my $tr1 ( grep { $_->[0] !~ m!^(chr|chrID|chromID)$!i } &fileSunhh::load_tabFile($opts{'chr_keyLen'}) ) {
		defined $chr_start{$tr1->[0]} and next; 
		$chr_start{$tr1->[0]} = $prevE+1; 
		$prevE = $prevE + $tr1->[1] + $opts{'insertLen'}; 
	}
}


#my $ofh_geno = &openFH( "$opts{'opref'}.geno", '>' ); 
#my $ofh_snp  = &openFH( "$opts{'opref'}.snp",  '>' ); 
#my $ofh_log  = &openFH( "$opts{'opref'}.pos_log", '>' ); 
my $para_colMax; 
my $in_cnt = 0; 
my $in_wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
my @in_toMerge; 
for my $fh (@InFp) {
	$in_cnt ++; 
	&tsmsg("[Msg] Processing [$in_cnt] [$fh]\n"); 
	# Deal with header_txt : Nothing to do. 
	# Sep-files
	my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
	my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" ); 
	# Process-sub-files
	### Setup %chr_start
	if ( !(defined $opts{'chr_keyLen'}) ) {
		for my $sfn ( @sub_fn ) {
			my $pid = $pm->start and next; 
			open F,'<',"$sfn" or die; 
			my %max_p; 
			my %ord; 
			while (<F>) {
				m/^\s*(#|$)/ and next;
				chomp; 
				my @ta = split(/\t/, $_); 
				if ($. == 1 and $ta[0] =~ m/^(chr|chrID|chromID)$/i) {
					next; 
				}
				$max_p{$ta[0]} //= $ta[1]; 
				$max_p{$ta[0]} < $ta[1] and $max_p{$ta[0]} = $ta[1]; 
				$ord{$ta[0]} //= $.; 
			}
			close F; 
			open O,'>',"$sfn.max_p" or die; 
			for my $tk (sort {$ord{$a} <=> $ord{$b}} keys %ord) {
				print O join("\t", $tk, $max_p{$tk})."\n"; 
			}
			close O; 
			$pm->finish; 
		}
		$pm->wait_all_children; 
		my %max_p; 
		my %ord; 
		my $tcnt = 0; 
		for my $sfn ( @sub_fn ) {
			open F,'<',"$sfn.max_p" or die; 
			while (<F>) {
				$tcnt ++; 
				chomp; 
				my @ta = split(/\t/, $_); 
				$ord{$ta[0]} //= $tcnt; 
				$max_p{$ta[0]} //= $ta[1]; 
				$max_p{$ta[0]} < $ta[1] and $max_p{$ta[0]} = $ta[1]; 
			}
			close F; 
		}
		my $prevE = 0; 
		for my $tk (sort {$ord{$a} <=> $ord{$b}} keys %ord) {
			defined $chr_start{$tk} and die "Strange [$tk]\n"; 
			$chr_start{$tk} = $prevE+1; 
			$prevE = $prevE + $max_p{$tk} + $opts{'insertLen'}; 
		}
	}
	### Generate output files. 
	for my $sfn ( @sub_fn ) {
		my $pid = $pm->start and next; 
		my ($ref_base, $var_base); 
		open F,'<',"$sfn" or die; 
		open O_SNP,'>',"$sfn.snp" or die; 
		open O_LOG,'>',"$sfn.pos_log" or die; 
		open O_GEN,'>',"$sfn.geno" or die; 
		LINE: 
		while (<F>) {
			m/^\s*(#|$)/ and next LINE; 
			chomp; 
			my @ta = split(/\t/, $_); 
			$para_colMax //= $#ta; 
			if ($. == 1 and $ta[0] =~ m/^(chr|chrID|chromID)$/i) { 
				next LINE; 
			}

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
			print O_LOG join("\t", "$ta[0]_MRK_$ta[1]", $ta[0], $ta[1], 1, $ta[1]+$chr_start{$ta[0]}-1)."\n"; 
			print O_SNP join("\t", "$ta[0]_MRK_$ta[1]", 1, '0.0', ($ta[1]+$chr_start{$ta[0]}-1), $ref_base, $var_base)."\n"; 
			
			for (my $i=$opts{'startColN'}; $i<=$para_colMax; $i++) {
				if ( $opts{'input_is_tab'} ) {
					$ta[$i] eq './.' and do { print O_GEN '9'; next; }; 
				} else {
					$ta[$i] eq 'N'   and do { print O_GEN '9'; next; }; 
				}
				my $cnt = 0; 
				for my $tb ( @{$tmp_al[$i]}[0,1] ) {
					$tb eq $ref_base and $cnt++; 
				}
				print O_GEN $cnt;
			}
			print O_GEN "\n"; 
		}
		close O_SNP; 
		close O_GEN; 
		close O_LOG; 
		close F; 
		$pm->finish; 
	}
	$pm->wait_all_children; 
	# Merge sub-files
	&merge_output( \@sub_fn, "$in_wrk_dir/$in_cnt" ); 
	push(@in_toMerge, "$in_wrk_dir/$in_cnt"); 
	# Delete temp_dir
	&fileSunhh::_rmtree($wrk_dir); 
}
if ($in_cnt >= 1) {
	&merge_output( \@in_toMerge, "$opts{'opref'}" ); 
} else {
	die "why?\n"; 
}
&fileSunhh::_rmtree($in_wrk_dir); 


sub merge_output {
	my $ipref_aref = shift; 
	my $opref = shift; 
	defined $opref or die "$_\n"; 
	if (@$ipref_aref == 1) {
		my $sfn = $ipref_aref->[0]; 
		&fileSunhh::_copy( "$sfn.snp", "$opref.snp" ); 
		&fileSunhh::_copy( "$sfn.geno", "$opref.geno" ); 
		&fileSunhh::_copy( "$sfn.pos_log", "$opref.pos_log" ); 
		return 0; 
	}

	open MO_SNP,'>',"${opref}.snp" or die; 
	open MO_GEN,'>',"${opref}.geno" or die; 
	open MO_LOG,'>',"${opref}.pos_log" or die; 
	my @sub_fn = @$ipref_aref; 
	my @geno_lines; 
	for my $sfn ( @sub_fn ) {
		open F_SNP,'<',"$sfn.snp" or die; 
		while (<F_SNP>) {
			print MO_SNP $_; 
		}
		close F_SNP; 
		open F_LOG,'<',"$sfn.pos_log" or die; 
		while (<F_LOG>) {
			print MO_LOG $_; 
		}
		close F_LOG; 
		open F_GEN,'<',"$sfn.geno" or die; 
		while (<F_GEN>) {
			print MO_GEN $_; 
		}
		close F_GEN; 
	}
	close MO_SNP; 
	close MO_GEN; 
	close MO_LOG; 
	return 0; 
}# merge_output () 
