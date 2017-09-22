#!/usr/bin/perl
# Reference : http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed 
# The autosomes should be coded 1 through 22. The following other codes can be used to specify other chromosome types:
#  X    X chromosome                    -> 23
#  Y    Y chromosome                    -> 24
#  XY   Pseudo-autosomal region of X    -> 25
#  MT   Mitochondrial                   -> 26
#
# All markers should be biallelic 
# Add more restrictions for genotypes: 
#  1. Only accept A/T/G/C or bi-allelic genotypes. 
#      Allow '*' as deletion , '+' as insertion. 
#      Do not allow '.' or '-' in sMao.cols format. 
#      Allow .vcf.tab format input. 
#  2. Only bi-allelic site accepted when using -max2allele. 
# 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	'help!', 
	'max2allele!', 
	'startColN:i', 
	'noHeader!', 
	'file_list:s', 
	'list_h2famID:s', 
	'list_h2indID:s', 
	'list_chrID2num:s', 
); 
$opts{'startColN'} //= 2; 

sub usage {
	print STDERR <<HH; 

perl $0 in.cols
# Will write to files : in.cols.ped in.cols.map in.cols.info

-file_list      [''] If given, will replace the ARGV with files in file_list. 

-help 
-max2allele     Only accept sites with exact TWO alleles if given; Should be given for haploview. 
                  This is also required by plink's .ped file format. 
-startColN      [2] 0-indexed number of the first column for sample base. 
-noHeader       Given if input .cols file has no header. 

-list_h2famID   [filename] Optional. Format : header_ID \\t family_ID
-list_h2indID   [filename] Optional. Format : header_ID \\t indivdual_ID
-list_chrID2num [filename] Optional. Recommended. Format : chr_ID \\t number (101 .. 999999)


HH
	exit 1; 
}

!@ARGV and !(defined $opts{'file_list'}) and &usage(); 
$opts{'help'} and &usage(); 

if (defined $opts{'file_list'}) {
	@ARGV = (); 
	my $l_fh = &openFH($opts{'file_list'}, '<'); 
	while (<$l_fh>) {
		m/^\s*$/ and next; 
		m/^\s*#/ and next; 
		chomp; 
		my @ta = split(/\t/, $_); 
		push(@ARGV, $ta[0]); 
	}
	close($l_fh); 
}

my %dblist; 
{
my @aa = (
[qw/W A T/], 
[qw/S C G/],
[qw/M A C/],
[qw/K G T/],
[qw/R A G/],
[qw/Y C T/]
); 
for my $tr (@aa) {
	my @bb = @$tr; 
	$dblist{$bb[0]} = [$bb[1], $bb[2]]; 
}
}

# .ped file format : 
#   pedID indvID FatherID MotherID sexTag caseTag genotypes(A=1, C=2, G=3, T=4, Miss=0)
my %allele2num = qw(
A 1
C 2
G 3
T 4
N 0
); 
my %num2allele = qw(
1 A
2 C
3 G
4 T
0 N
); 
my %glob; 
$glob{'al2num_max'} = (sort { $b <=> $a } keys %num2allele)[0]; 

my (%h2fam, %h2ind, %used, %chr2num, %num2chr); 
if ( defined $opts{'list_h2famID'} ) {
	%h2fam = map { $_->[0] => $_->[1] } &fileSunhh::load_tabFile( $opts{'list_h2famID'} ); 
}
if ( defined $opts{'list_h2indID'} ) {
	%h2ind = map { $_->[0] => $_->[1] } &fileSunhh::load_tabFile( $opts{'list_h2indID'} ); 
}
if ( defined $opts{'list_chrID2num'} ) {
	%chr2num = map { $_->[0] => $_->[1] } &fileSunhh::load_tabFile( $opts{'list_chrID2num'} ); 
	for my $tk (keys %chr2num) {
		defined $num2chr{$chr2num{$tk}} and &stopErr("[Err] Repeat chr-number [$chr2num{$tk}]\n"); 
		$num2chr{ $chr2num{$tk} } = $tk; 
	}
}



for my $inF (@ARGV) {

	&tsmsg("[Rec] Processing [$inF]\n"); 

	my $oPedFile = "$inF.ped"; 
	my $oMapFile = "$inF.map"; 
	my $oInfFile = "$inF.info"; 
	
	my $in_fh = &openFH($inF, '<'); 
	open OM,'>',"$oMapFile" or die; 
	open OI,'>',"$oInfFile" or die; 
	# .map format : chrID, SNP_ID, cM_num(0), position
	my (@header, @data); 
	my %tmp_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>5e6 ); 
	SNP_LINE: 
	while (<$in_fh>) {
		&fileSunhh::log_section( $. , \%tmp_cnt ) and &tsmsg("[Msg]   Loading $. line in file [$inF].\n"); 
		chomp; 
		my @ta = split(/\t/, $_); 
		if ($. == 1 and !($opts{'noHeader'})) {
			@header = @ta; 
			next; 
		}
		my %cnt_allele; 
		my @tmp_nn; 
		for (my $i=$opts{'startColN'}; $i<@ta; $i++) {
			my @nums = &geno2num($ta[$i]); 
			if ($nums[0] eq '') {
				### Check if this site is good. 
				defined $used{'bad_geno'}{$ta[$i]} or &tsmsg("[Wrn] Skip line with genotype [$ta[$i]]\n"); 
				$used{'bad_geno'}{$ta[$i]} = 1; 
				next SNP_LINE; 
			}
			for my $tn (@nums) {
				$tn == 0 and next; 
				$cnt_allele{$tn} ++; 
			}
			$tmp_nn[$i] = [@nums]; 
		}
		if ($opts{'max2allele'}) {
			scalar(keys %cnt_allele) == 2 or next SNP_LINE; 
		} else {
			scalar(keys %cnt_allele) <= 1 and next SNP_LINE; 
		}
		my @srt_allele = sort { $cnt_allele{$b} <=> $cnt_allele{$a} || $a cmp $b } keys %cnt_allele; 
		my %use; 
		for my $tn (@srt_allele[0,1]) {
			$use{$tn} = 1; 
		}
		for ( my $i=$opts{'startColN'}; $i<@tmp_nn; $i++) {
			my $is_0 = 0; 
			for my $tn (@{$tmp_nn[$i]}) {
				defined $use{$tn} or do { $is_0 = 1; last; }; 
			}
			$is_0 == 1 and @{$tmp_nn[$i]} = (0,0); 
			push(@{$data[$i]}, "$tmp_nn[$i][0] $tmp_nn[$i][1]"); 
		}
		my $chrID = $ta[0]; 
		if ( defined $chr2num{$chrID} ) {
			$chrID = $chr2num{$chrID}; 
		} else {
			my $raw_chrID = $chrID; 
			$chrID =~ s/^chr(\d+)$/$1/i; $chrID =~ s/^0+$/0/; 
			$chrID =~ s/^WM97_Chr0*(\d*)$/$1/i; 
			$chrID eq '' and $chrID = 999999; 
			defined $num2chr{$chrID} and &stopErr("[Err] Repeat using chr-number [$chrID] for $raw_chrID\n"); 
			$chr2num{$raw_chrID} = $chrID; 
			$num2chr{$chrID} = $raw_chrID; 
		}
		print OM join("\t", $chrID, "s${chrID}_$ta[1]", 0, $ta[1])."\n"; 
		my $af_0 = sprintf("%.4f", $cnt_allele{$srt_allele[0]}/( $cnt_allele{$srt_allele[0]} + $cnt_allele{$srt_allele[1]}) ); 
		my $af_1 = sprintf("%.4f", $cnt_allele{$srt_allele[1]}/( $cnt_allele{$srt_allele[0]} + $cnt_allele{$srt_allele[1]}) ); 
		print OI join("\t", $#{$data[$opts{'startColN'}]}+1, $ta[1], $num2allele{ $srt_allele[0] }, $af_0, $num2allele{ $srt_allele[1] }, $af_1)."\n"; 
	}
	close($in_fh); 
	close OM; 
	close OI; 
	# output .ped 
	open OP,'>',"$oPedFile" or die; 
	for (my $i=0; $i+$opts{'startColN'}<@data; $i++) {
		my $pedID = ($opts{'noHeader'}) ? $i+1 : $header[$i+$opts{'startColN'}]; 
		my $idvID = ($opts{'noHeader'}) ? $i+1 : $header[$i+$opts{'startColN'}]; 
		defined $h2fam{$pedID} and $pedID = $h2fam{$pedID}; 
		defined $h2ind{$idvID} and $idvID = $h2ind{$idvID}; 
		my $comb_ID = "${pedID}\t$idvID"; 
		defined $used{'comb_ID'}{$comb_ID} and &stopErr("[Err] Repeated ID [$comb_ID]\n"); 
		$used{'comb_ID'}{$comb_ID} = 1; 
		my $fatID = 0; 
		my $monID = 0; 
		my $sexID = 0; # 0/1? I used 1 at first, but found BGI use 0 instead. 0=male, 2=female, other=unknown; 
		my $case  = 0; # -9 missing ; 0 missing ; 1 unaffected ; 2 affected. 
		my $genotypeLine = join("\t", @{$data[$i+$opts{'startColN'}]}); 
		print OP join("\t", $pedID, $idvID, $fatID, $monID, $sexID, $case, $genotypeLine)."\n"; 
	}
	close OP; 
	
}
&tsmsg("[Rec] All done.\n"); 

sub get_al2num {
	unless ( defined $allele2num{$_[0]} ) {
		$glob{'al2num_max'} ++; 
		$allele2num{ $_[0] } = $glob{'al2num_max'}; 
		$num2allele{ $glob{'al2num_max'} } = $_[0]; 
	}
	return( $allele2num{$_[0]} ); 
}# get_al2num() 

sub geno2num {
	my $tb = shift; 
	if ( $tb =~ m!^([^/\s]+)/([^/\s]+)$!i ) {
		# .vcf.tab format ; 
		my ($a1, $a2) = (uc($1), uc($2)); 
		if ( $a1 eq '.' or $a1 eq 'N' ) {
			$a1 eq $a2 or &stopErr("[Err] Bad genotype format [$tb]\n"); 
			return( &get_al2num( 'N' ) , &get_al2num( 'N' ) ); 
		} elsif ( $a1 =~ m/^[ATGCN*]$/i and $a2 =~ m/^[ATGCN*]$/i ) {
			my $n1 = &get_al2num($a1) ; 
			my $n2 = &get_al2num($a2) ; 
			return( $n1 , $n2 ); 
		} else {
			my $n1 = &get_al2num( $a1 ); 
			my $n2 = &get_al2num( $a2 ); 
			return( $n1, $n2 ); 
		}
		# ( $a1 =~ m/^[ATGCN*]$/i and $a2 =~ m/^[ATGCN*]$/i ) or &stopErr("[Err] genotype_2 [$tb]\n"); 
	}
	# Mao's pipeline output SNP table format: 
	$tb =~ m/^N+$/i and return( &get_al2num( 'N' ) , &get_al2num( 'N' ) ); 
	if ( $tb =~ m/^[ATGC\*]$/i ) {
		my $n1 = &get_al2num( uc($tb) ); 
		return( $n1, $n1 ); 
	} elsif ( $tb =~ m/^([ATGC\*])([ATGC\*])$/i ) {
		my ($a1, $a2) = (uc($1), uc($2)); 
		my $n1 = &get_al2num( $a1 ); 
		my $n2 = &get_al2num( $a2 ); 
		return( $n1, $n2 ); 
	} elsif ( defined $dblist{ uc($tb) } ) {
		my @b1 = @{ $dblist{ uc($tb) } }; 
		my $n1 = &get_al2num( uc($b1[0]) ); 
		my $n2 = &get_al2num( uc($b1[1]) ); 
		return( $n1,$n2 ); 
	} elsif ( $tb =~ m/^[ATGCN]\+([ATGCN]+)$/i ) {
		my $n1 = &get_al2num( uc("+$1") ); 
		return( $n1, $n1 ); 
	} elsif ( $tb =~ m/^([ATGCN])[ATGCN]\+([ATGCN]+)$/i ) {
		my $a1 = uc($1); 
		my $a2 = uc("+$2"); 
		my $n1 = &get_al2num( $a1 ); 
		my $n2 = &get_al2num( $a2 ); 
		return( $n1,$n2 ); 
	} else {
		return( '', '' ); 
	}
}

