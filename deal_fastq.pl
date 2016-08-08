#!/usr/bin/env perl
# 2014-03-20 A script to deal with fastq format reads. It will be always in processing. 
# 2014-03-25 Add function to search a pattern in reads, and return different values as defined. 
#

BEGIN {
	use File::Basename; 
	use Cwd 'abs_path'; 
	my $pl_dir = dirname( abs_path(__FILE__) ); 
	use lib "$pl_dir/MyPM"; 
}


use strict; 
use warnings; 
use List::Util qw(first max maxstr min minstr reduce shuffle sum); 
use Getopt::Long; 
use LogInforSunhh; 
use fileSunhh; 
my %opts; 

sub usage {
	print STDOUT <<HELP; 
#******* Instruction of this program *********#

Introduction:Deal with fasta format file. 
perl $0 in.fastq 

-fq2fa        [Boolean] Transformat fastq sequences to fasta format. 
-oQfile       [String] Name of quality file to output for -fq2fa fastq sequence. 

-fq2Val       [Boolean] Change the quality string to a set of quality numbers seprated by white space. 

-keep_len     [String] "min_len-max_len". Extract sequences whose lengths are between min_len and max_len.
-paired       [Boolean] Treat files as paired, and output reads in order of R1/R2 ... 

-sepR12_pref  [String] Input reads as R1-R2 order, and we output R1 to pref_R1.fq and R2 to pref_R2.fq
-joinR12      [Boolean] Input files as paired, and join R1/R2 reads in a same stream. You can use \"-keep_len 0- -paired\" to get the same result. 

-showQscale   [Boolean] Show fastq quality scale value [33/64]. 

-33to64       [Boolean] Transform phred33 to phred64 quality. 
-64to33       [Boolean] Transform phred64 to phred33 quality. 

-rd_Num       [Boolean] Summary read number/length of fq files. 
-rd_Num_fast  [Boolean] Read simple format of fastq faster. Need -rd_Num assigned. 

-rd_LenHist   [Boolean] 
-rd_LenHist_pmin   [Boolean] by paired-min. 
-rd_LenHist_range  [1/0-max] Example: [10-20], then 10 stands for <=10, 20 stands for >=20; 
-rd_LenHist_name   ["RdNum"] A string for the 2nd column's name. 

-frag         [String] [start-end] output a fragment of the squence in single sequence fasta file;Start position is 1; May be s1-e1:s2-e2...
-frag_r       [Boolean] give out reversed string; 
-frag_c       [Boolean] give out complemented string as DNA; 

-phred_scale  [INT] 33/64. Set scale value to get phred quality score. 
-guess_scale  [Boolean] Guess the phred scale if given. This will overwrite -phred_scale option. 

-search       [String] Pattern to search in read. 
-srch_strand  [f/r/b] forward/reverse/both. 
-srch_back    [position/match/both/read] Default read
-srch_drop    [Boolean] Only valid if -srch_back == read. It will drop the reads matching to the pattern. 
-srch_max     [Boolean] Return all patterns bp by bp if given. 

-get_key      [Boolean] Return read keys



-rdKey        [Boolean] Keep only the first non-blank characters for read key if given. 

-randSlct     [1.0] Random select reads/pairs to a subset of ratio from (0-1]. 

-reorder      [orderList] Format : rdID_01 \\n rdID_02 \\n rdID_03 \\n ...

#******* Instruction of this program *********#
HELP
	exit(1); 
}

GetOptions(\%opts, 
	"phred_scale:i", "guess_scale!", 
	"rdKey!", 
	"fq2fa!", "oQfile:s", 
	"fq2Val!", 
	"keep_len:s", "paired!", 
	"sepR12_pref:s", 
	"joinR12!", # In fact, this function can be performed by keep_len() subroutine with "-keep_len 0- -paired "
	"showQscale!", 
	"33to64!", "64to33!", 
	"rd_Num!", "rd_Num_fast!", "rd_LenHist!", "rd_LenHist_pmin!", "rd_LenHist_range:s", "rd_LenHist_name:s", 
	"frag:s", "frag_r!", "frag_c!", 
	"search:s", "srch_strand:s", "srch_back:s", "srch_drop!", "srch_max!", 
	"randSlct:f", 
	"get_key!", 
	"reorder:s", 
	"help!", 
); 

#****************************************************************#
#--------------Main-----Function-----Start-----------------------#
#****************************************************************#

&usage() if ( $opts{help} or ( -t and !@ARGV ) );
$opts{'phred_scale'} //= 33; 
defined $opts{'guess_scale'} and $opts{'phred_scale'} = undef(); 


# Making File handles for reading;
our @InFp = () ; # 2007-8-29 16:07 È«¾Ö±äÁ¿!
if ( !@ARGV )
{
	@InFp = (\*STDIN);
}
else
{
	for (@ARGV) {
		push( @InFp, &openFH($_,'<') );
	}
}


## Global settings
my %good_str = qw(
	f        1
	r        -1
	b        2
	foward 	 1
	reverse  -1
  both     2
  +        1
  -        -1
  1        1
  -1       -1
  2        2
); 


## Call functions. 
&fq2fa() if ( $opts{fq2fa} ); 
&keep_len( $opts{keep_len} ) if ( defined $opts{keep_len} and $opts{keep_len} ne '' ); 
&showQscale() if ( $opts{showQscale} ); 
&fq2fq( 64-33 ) if ( $opts{'33to64'} ); 
&fq2fq( 33-64 ) if ( $opts{'64to33'} ); 
&sepR12( $opts{sepR12_pref} ) if ( defined $opts{sepR12_pref} and $opts{sepR12_pref} ne '' ); 
&joinR12() if ( $opts{joinR12} ); 
&fq2Val() if ( $opts{fq2Val} ); 
&rd_Num() if ( $opts{rd_Num} ); 
&fragmentRd() if ( defined $opts{frag} and $opts{frag} ne '' ); 
&searchPattern() if ( defined $opts{search} and $opts{search} ne '' ); 
&rd_LenHist() if ( $opts{rd_LenHist} ); 
&randSlct() if ( defined $opts{'randSlct'} ); 
&get_key() if ( $opts{'get_key'} ); 
&reorder() if ( defined $opts{'reorder'} ); 

#****************************************************************#
#--------------Subprogram------------Start-----------------------#
#****************************************************************#

# -randSlct     [1]
# -paired       
sub randSlct {
	$opts{'randSlct'} >= 1 and &stopErr("[Err] No need to select.\n"); 
	$opts{'randSlct'} <= 0 and &stopErr("[Err] Bad setting [$opts{'randSlct'}]\n"); 
	if ( $opts{'paired'} ) {
		for (my $i=0; $i<@InFp; $i+=2) {
			my $fh1 = $InFp[$i];
			my $fh2 = $InFp[$i+1];
			RD: 
			while ( !eof($fh1) and !eof($fh2) ) {
				my $rdRec1 = &get_fq_record($fh1);
				my $rdRec2 = &get_fq_record($fh2);
				rand(1) <= $opts{'randSlct'} or next; 
				print STDOUT "\@$rdRec1->{id}$rdRec1->{seq}+\n$rdRec1->{qual}\n"; 
				print STDOUT "\@$rdRec2->{id}$rdRec2->{seq}+\n$rdRec2->{qual}\n"; 
			}#End while() RD:
		}
	} else {
		for (my $i=0; $i<@InFp; $i++) {
			my $fh1 = $InFp[$i]; 
			RD: 
			while ( !eof($fh1) ) {
				my $rdRec1 = &get_fq_record($fh1);
				rand(1) <= $opts{'randSlct'} or next; 
				print STDOUT "\@$rdRec1->{id}$rdRec1->{seq}+\n$rdRec1->{qual}\n"; 
			}#End while() RD:
		}
	}
}#End randSlct() 

# -search       [String] Pattern to search in read. 
# -srch_strand  [f/r/b] forward/reverse/both. 
# -srch_back    [position/match/both/read] Default read
# -srch_drop    [Boolean] Only valid if -srch_back == read. It will drop the reads matching to the pattern. 
sub searchPattern {
	my $srch_drop = 0; defined $opts{srch_drop} and $srch_drop = 1; 
	my $srch_str = $good_str{b}; ; defined $opts{srch_strand} and $srch_str = $good_str{ $opts{srch_strand} }; 
	defined $good_str{ $srch_str } or &stopErr( "[Err] Failed to get strand information! -srch_strand == [$srch_str]\n" ); 
	my $srch_back = 'read'; defined $opts{srch_back} and $srch_back = $opts{srch_back}; 
	my $srch_mm = 'min'; defined $opts{srch_max} and $srch_mm = 'max'; 
	my $srch_pattern = $opts{search}; 
	if ( $srch_back eq 'position' ) {
		print STDOUT join("\t", qw/rdID rdLen MatchStart MatchEnd MatchLen/)."\n"; 
	} elsif ( $srch_back eq 'match' ) {
		print STDOUT join("\t", qw/rdID rdLen MatchLen MatchSeq/)."\n"; 
	} elsif ( $srch_back eq 'both' ) {
		print STDOUT join("\t", qw/rdID rdLen MatchStart MatchEnd MatchLen MatchSeq/)."\n"; 
	} elsif ( $srch_back eq 'read' ) {
		; 
	} else {
		&stopErr("[Err] Unkonwn -srch_back [$srch_back]\n"); 
	}

	if ( $srch_str == 1 ) {
		for my $fh ( @InFp ) {
			while ( my $rdRec = &get_fq_record($fh) ) {
				$rdRec->{seq} =~ s/\s//g; 
				my $rdLen = length($rdRec->{seq}); 
				my @t_match = &siteList( \$srch_pattern, \$rdRec->{seq}, $srch_mm ); 
				if ( $srch_back eq 'read' ) {
					if ( scalar(@t_match) > 0 ) {
						$srch_drop == 0 and print STDOUT "\@$rdRec->{id}$rdRec->{seq}\n+\n$rdRec->{qual}\n"; 
					} else {
						$srch_drop == 1 and  print STDOUT "\@$rdRec->{id}$rdRec->{seq}\n+\n$rdRec->{qual}\n"; 
					}
				} elsif ( $srch_back eq 'position' ) {
					chomp($rdRec->{id}); 
					map { print STDOUT join("\t", $rdRec->{id}, $rdLen, $_->[0], $_->[1], $_->[1]-$_->[0]+1)."\n"; } @t_match; 
				} elsif ( $srch_back eq 'match' ) {
					chomp($rdRec->{id}); 
					map { print STDOUT join("\t", $rdRec->{id}, $rdLen, $_->[1]-$_->[0]+1, $_->[2])."\n"; } @t_match; 
				} elsif ( $srch_back eq 'both' ) {
					chomp($rdRec->{id}); 
					map { print STDOUT join("\t", $rdRec->{id}, $rdLen, $_->[0], $_->[1], $_->[1]-$_->[0]+1, $_->[2])."\n"; } @t_match; 
				} else {
					&stopErr("[Err] Unknown -srch_back [$srch_back]\n"); 
				}
			}#End while ()
		}#End for my $fh
	} elsif ( $srch_str == -1 ) {
		for my $fh ( @InFp ) {
			while ( my $rdRec = &get_fq_record($fh) ) {
				$rdRec->{seq} =~ s/\s//g; 
				my $rdLen = length($rdRec->{seq}); 
				my $seqRC = $rdRec->{seq}; &rcSeq( \$seqRC, 'rc' ); 
				my @t_match = &siteList( \$srch_pattern, \$seqRC, $srch_mm ); 
				if ( $srch_back eq 'read' ) {
					if ( scalar(@t_match) > 0 ) {
						$srch_drop == 0 and print STDOUT "\@$rdRec->{id}$rdRec->{seq}\n+\n$rdRec->{qual}\n"; 
					} else {
						$srch_drop == 1 and print STDOUT "\@$rdRec->{id}$rdRec->{seq}\n+\n$rdRec->{qual}\n"; 
					}
				} elsif ( $srch_back eq 'position' ) {
					chomp($rdRec->{id}); 
					map { print STDOUT join("\t", $rdRec->{id}, $rdLen, $rdLen - $_->[0] + 1, $rdLen - $_->[1] + 1, $_->[1]-$_->[0]+1)."\n"; } @t_match; 
				} elsif ( $srch_back eq 'match' ) {
					chomp($rdRec->{id}); 
					map { &rcSeq(\$_->[2], 'rc'); print STDOUT join("\t", $rdRec->{id}, $rdLen, $_->[1]-$_->[0]+1, $_->[2])."\n"; } @t_match; 
				} elsif ( $srch_back eq 'both' ) {
					chomp($rdRec->{id}); 
					map { &rcSeq(\$_->[2], 'rc'); print STDOUT join("\t", $rdRec->{id}, $rdLen, $rdLen - $_->[0] + 1, $rdLen - $_->[1] + 1, $_->[1]-$_->[0]+1, $_->[2])."\n"; } @t_match; 
				} else {
					&stopErr("[Err] Unknown -srch_back [$srch_back]\n"); 
				}
			}#End while ()
		}#End for my $fh
	} elsif ( $srch_str == 2 ) {
		for my $fh ( @InFp ) {
			while ( my $rdRec = &get_fq_record($fh) ) {
				$rdRec->{seq} =~ s/\s//g; 
				my $seqRC = $rdRec->{seq}; &rcSeq( \$seqRC, 'rc' ); 
				my $rdLen = length($rdRec->{seq}); 
				my @t_match = &siteList( \$srch_pattern, \$rdRec->{seq}, $srch_mm ); 
				my @t_match_RC = &siteList( \$srch_pattern, \$seqRC, $srch_mm ); 
				if ( $srch_back eq 'read' ) {
					if ( scalar(@t_match) > 0 or scalar(@t_match_RC) > 0 ) {
						$srch_drop == 0 and print STDOUT "\@$rdRec->{id}$rdRec->{seq}\n+\n$rdRec->{qual}\n"; 
					} else {
						$srch_drop == 1 and print STDOUT "\@$rdRec->{id}$rdRec->{seq}\n+\n$rdRec->{qual}\n"; 
					}
				} elsif ( $srch_back eq 'position' ) {
					chomp($rdRec->{id}); 
					map { print STDOUT join("\t", $rdRec->{id}, $rdLen, $_->[0], $_->[1], $_->[1]-$_->[0]+1)."\n"; } @t_match; 
					map { print STDOUT join("\t", $rdRec->{id}, $rdLen, $rdLen - $_->[0] + 1, $rdLen - $_->[1] + 1, $_->[1]-$_->[0]+1)."\n"; } @t_match_RC; 
				} elsif ( $srch_back eq 'match' ) {
					chomp($rdRec->{id}); 
					map { print STDOUT join("\t", $rdRec->{id}, $rdLen, $_->[1]-$_->[0]+1, $_->[2])."\n"; } @t_match; 
					map { &rcSeq(\$_->[2], 'rc'); print STDOUT join("\t", $rdRec->{id}, $rdLen, $_->[1]-$_->[0]+1, $_->[2])."\n"; } @t_match_RC; 
				} elsif ( $srch_back eq 'both' ) {
					chomp($rdRec->{id}); 
					map { print STDOUT join("\t", $rdRec->{id}, $rdLen, $_->[0], $_->[1], $_->[1]-$_->[0]+1, $_->[2])."\n"; } @t_match; 
					map { &rcSeq(\$_->[2], 'rc'); print STDOUT join("\t", $rdRec->{id}, $rdLen, $rdLen - $_->[0] + 1, $rdLen - $_->[1] + 1, $_->[1]-$_->[0]+1, $_->[2])."\n"; } @t_match_RC; 
				} else {
					&stopErr("[Err] Unknown -srch_back [$srch_back]\n"); 
				}
			}#End while ()
		}#End for my $fh
	} else {
		&stopErr("[Err] Unknown -srch_strand [$srch_str]\n"); 
	}

}#sub searchPattern 

# Input : ( \$qry_pattern, \$sbjct_seq, 'min/max')
# Return: ( [s1,e1,matchSeq1], [s2,e2,matchSeq2], ... )
sub siteList ($$$) {
	my $siteR = shift or &stopErr( "[Err]no Site input!\n" );
	if ($$siteR eq "") {
		warn "[Err]Site sequence cannot be empty!\n";
		return;
	}
	my $qrSite = qr/$$siteR/s;
	my $refR = shift or &stopErr( "[Err]no refSeq input!\n" );
	my $modeChk = shift;
	my $Is_min = 1;
	if (defined $modeChk) {
		$modeChk = lc($modeChk);
		if ($modeChk eq 'max') {
			$Is_min = 0;
		} elsif ($modeChk eq 'min') {
			$Is_min = 1;
		} else {
			warn "[Err]Third para of \&siteList should be 'Max' or 'Min' . Not [$modeChk].\n";
		}
	}
	my @posList = ();
	pos($$refR) = 0;
	while ($$refR =~ m/\G(?:.*?)($qrSite)/gs) {
		push ( @posList, [ $-[1]+1 , $+[1], $1 ] );
		$Is_min?(pos($$refR) = $+[1]):(pos($$refR) = $-[1]+1);
	}
	return @posList;
}# end siteList subroutine. 2013-10-30



sub fragmentRd{
	my $fragR = &frag_regions($opts{frag}); 
	for my $fh ( @InFp ) {
		while ( my $rdRec = &get_fq_record($fh) ) {
			$rdRec->{seq} =~ s/\s//g; 
			$rdRec->{qual} =~ s/\s//g; 
			my $l_seq = length($rdRec->{seq}); 
			my ($strSeq, $strQual) = ('', ''); 
			my @Range; 
			for my $seR ( @$fragR ) {
				my ($add_s, $add_e) = @$seR; 
				$add_s < 0 and $add_s = $l_seq+$add_s+1; 
				$add_e eq 'end' and $add_e = $l_seq; 
				$add_e < 0 and $add_e = $l_seq+$add_e+1; 
				$add_e > $l_seq and $add_e = $l_seq; 
				if ( $add_s > $l_seq ) {
					; # Nothing happening. 
				} else {
					$strSeq  .= substr($rdRec->{seq},  $add_s-1, $add_e-$add_s+1); 
					$strQual .= substr($rdRec->{qual}, $add_s-1, $add_e-$add_s+1); 
				}
				push(@Range, "$add_s\-$add_e"); 
			}
			my $range = join(',', @Range); 
			if ( $opts{frag_c} ) { &rcSeq(\$strSeq, 'c'); $range = "C$range"; } 
			if ( $opts{frag_r} ) { $strSeq = reverse($strSeq); $strQual = reverse($strQual); $range = "R$range"; } 
			chomp($rdRec->{id}); $rdRec->{id} .= " [$range]"; 
			print STDOUT "\@$rdRec->{id}\n$strSeq\n+\n$strQual\n"; 
		}#while
	}#for
}#sub fragmentRd() 

sub showQscale {
	print STDOUT join("\t", qw/InFqName QualScale/)."\n"; 
	if ( !@ARGV ) {
		my $qscale = &guessPhredScale($InFp[0]); 
		print STDOUT join("\t", '[STDIN]', $qscale)."\n"; 
	} else {
		for (my $i=0; $i<@ARGV; $i++) {
			my $qscale = &guessPhredScale($InFp[$i]); 
			print STDOUT join("\t", $ARGV[$i], $qscale)."\n"; 
		}
	}
}#sub showQscale() 

sub fq2fa {
	my $oQfh = undef(); 
	if (defined $opts{oQfile} and $opts{oQfile} ne '') {
		$oQfh = &openFH($opts{oQfile}, '>'); 
	}
	my $rdNum = 0; 
	for my $fh ( @InFp ) {
		my $qscale = &guessPhredScale($fh); 
		if ( defined $oQfh ) {
			while ( my $rdRec = &get_fq_record($fh) ) {
				$rdNum ++; 
				$rdNum % 10e6 == 1 and &tsmsg("[Msg] Dealing $rdNum reads.\n"); 
				print STDOUT ">$rdRec->{id}$rdRec->{seq}"; 
				my $qVR = &qChar2qValue(\$rdRec->{qual}, $qscale); 
				print {$oQfh} ">$rdRec->{id}@$qVR\n"; 
			}
		} else {
			while ( my $rdRec = &get_fq_record($fh) ) {
				$rdNum ++; 
				$rdNum % 10e6 == 1 and &tsmsg("[Msg] Dealing $rdNum reads.\n"); 
				print STDOUT ">$rdRec->{id}$rdRec->{seq}"; 
			}
		}
	}# End for my $fh
	&tsmsg("[Rec] Finish $rdNum reads.\n"); 
}#sub fq2fa

sub reorder {
	my %rd_rec; 
	my %cnt = ( 'cur_ln' => 0 , 'cntN_step' => 5e6 ); 
	for my $fh ( @InFp ) {
		while ( my $rdRec = &get_fq_record($fh) ) {
			$cnt{'cur_ln'} ++; 
			&fileSunhh::log_section( $cnt{'cur_ln'}, \%cnt ) and &tsmsg("[Msg] Dealing with $cnt{'cur_ln'} reads.\n"); 
			chomp( $rdRec->{'id'} ); 
			$rd_rec{ $rdRec->{'id'} } //= "\@$rdRec->{'ID'}$rdRec->{'seq'}+\n$rdRec->{'qual'}\n"; 
		}
	}

	my $fh_order = &openFH( $opts{'reorder'}, '<' ); 
	%cnt = ( 'cur_ln' => 0, 'cntN_step' => 5e6 ); 
	$cnt{'errCnt'} = 0; 
	while (<$fh_order>) {
		&fileSunhh::log_section( $. , \%cnt ) and &tsmsg("[Msg] Dealing with $cnt{'cur_ln'} line.\n"); 
		chomp; 
		my @ta=split(/\t/, $_); 
		unless ( defined $rd_rec{ $ta[0] } ) {
			$cnt{'errCnt'} ++; 
			$cnt{'errCnt'} > 50 or &tsmsg("[Wrn] Skip missing record [$ta[0]]\n"); 
			next; 
		}
		print STDOUT $rd_rec{ $ta[0] }; 
	}
	close($fh_order); 
	
	return; 
}# sub reorder () 

sub get_key {
	my %cnt; 
	$cnt{'cur_ln'} = 0; 
	$cnt{'cntN_step'} = 5e6; 
	for my $fh ( @InFp ) {
		while ( my $rdRec = &get_fq_record($fh) ) {
			$cnt{'cur_ln'} ++; 
			&fileSunhh::log_section( $cnt{'cur_ln'}, \%cnt ) and &tsmsg("[Msg] Dealing $cnt{'cur_ln'} reads.\n"); 
			print STDOUT $rdRec->{'id'}; 
		}
	}
	&tsmsg("[Rec] get_key() done for $cnt{'cur_ln'} reads.\n"); 
}# sub get_key ()

sub rd_Num {
	print STDOUT join("\t", qw/InFile Total_size Total_Rd_num Mean_Rd_size Range_Rd_size PhredCut Time/)."\n"; 
	for (my $i=0; $i<@InFp; $i++) {
		my $fh = $InFp[$i]; 
		my $fn = ( defined $ARGV[$i] ) ? $ARGV[$i] : '[STDIN]' ; 
		my $qscale = &guessPhredScale($fh); 
		&tsmsg("[Rec] Reading file [$fn] with phred scale [$qscale]\n"); 
		my %dd; 
		$dd{rdN} = $dd{bpN} = 0; 
		$dd{len_range} = ['NA', 'NA']; 
		unless ( $opts{'rd_Num_fast'} ) {
			while ( my $rdRec = &get_fq_record($fh) ) {
				$dd{rdN} ++; $dd{rdN} % 10e6 == 1 and &tsmsg("[Msg] $dd{rdN} reads.\n"); 
				$rdRec->{seq} =~ s/\s//g; 
				my $tt_len = length( $rdRec->{seq} ); 
				$dd{bpN} += $tt_len; 
				if ($dd{len_range}[0] eq 'NA') {
					$dd{len_range}[0] = $tt_len; 
					$dd{len_range}[1] = $tt_len; 
				} else {
					if      ( $dd{len_range}[0] > $tt_len ) {
						$dd{len_range}[0] = $tt_len; 
					} elsif ( $dd{len_range}[1] < $tt_len ) {
						$dd{len_range}[1] = $tt_len; 
					}
				}
			}#while
		} else {
			while ( <$fh> ) {
				$_ = <$fh>; chomp($_); <$fh>; <$fh>; 
				$dd{'rdN'} ++; $dd{'rdN'} % 10e6 == 1 and &tsmsg("[Msg] $dd{rdN} reads.\n"); 
				my $tt_len = length( $_ ); 
				$dd{'bpN'} += $tt_len; 
				if ( $dd{'len_range'}[0] eq 'NA' ) {
					$dd{'len_range'}[0] = $tt_len; 
					$dd{'len_range'}[1] = $tt_len; 
				} else {
					if ( $dd{'len_range'}[0] > $tt_len ) {
						$dd{'len_range'}[0] = $tt_len; 
					} elsif ( $dd{'len_range'}[1] < $tt_len ) {
						$dd{'len_range'}[1] = $tt_len; 
					}
				}
			}
		}
		$dd{avgL} = ( $dd{rdN} == 0 ) ? 'NA' : $dd{bpN}/$dd{rdN} ; 
		print STDOUT join("\t", $fn, $dd{bpN}, $dd{rdN}, $dd{avgL}, "$dd{len_range}[0]-$dd{len_range}[1]", "Phred$qscale", scalar(localtime()))."\n"; 
		&tsmsg("[Rec] Finish $dd{rdN} reads in file [$fn].\n"); 
	}#for 
}#sub rd_Num()



sub fq2Val {
	local $" = " "; 
	for my $fh ( @InFp ) {
		my $qscale = &guessPhredScale($fh); 
		while ( my $rdRec = &get_fq_record($fh) ) {
			my $qVR = &qChar2qValue(\$rdRec->{qual}, $qscale); 
			print STDOUT "\@$rdRec->{id}$rdRec->{seq}+\n@$qVR\n"; 
		}
	}
}#sub fq2Val() 

# -rd_LenHist
# -rd_LenHist_pmin
# -rd_LenHist_range : [1-maxLen]
# -rd_LenHist_name  : [rdNum]
sub rd_LenHist {
	my $is_pmin = 0; 
	my %lhist; 
	defined $opts{rd_LenHist_name} or $opts{rd_LenHist_name} = 'RdNum'; 

	if ( $opts{rd_LenHist_pmin} ) {
		for ( my $i=0; $i<@InFp; $i+=2 ) {
			my $fh1 = $InFp[$i]; 
			my $fh2 = $InFp[$i+1]; 
			while ( !eof($fh1) and !eof($fh2) ) {
				my $rdRec1 = &get_fq_record($fh1); 
				my $rdRec2 = &get_fq_record($fh2); 
				(my $ss1 = $rdRec1->{seq}) =~ s/\s//g; 
				(my $ss2 = $rdRec2->{seq}) =~ s/\s//g; 
				my $l1 = length($ss1); 
				my $l2 = length($ss2); 
				my $lmin = ( $l1 < $l2 ) ? $l1 : $l2 ; 
				$lhist{$lmin} ++; 
			}
		}
	} else {
		for my $fh ( @InFp ) {
			while ( my $rdRec = &get_fq_record($fh) ) {
				(my $ss = $rdRec->{seq}) =~ s/\s//g; 
				my $ll = length($ss); 
				$lhist{$ll} ++; 
			}
		}
	}
	my @srt_len = sort { $a <=> $b } keys %lhist; 
	my ($min, $max) = @srt_len[0, $#srt_len]; 
	$min < 1 or $min = 1; 
	if ( defined $opts{rd_LenHist_range} ) {
		$opts{rd_LenHist_range} =~ m/^(\d+)\-(\d+)$/ or &stopErr("[Err] -rd_LenHist_range unparsed [$opts{rd_LenHist_range}]\n"); 
		($min, $max) = ($1, $2); 
		my %new_lhist; 
		for my $ll ( keys %lhist ) {
			if ( $ll < $min ) {
				$new_lhist{$min} += $lhist{$ll}; 
			} elsif ( $ll > $max ) {
				$new_lhist{$max} += $lhist{$ll}; 
			} else {
				$new_lhist{$ll} += $lhist{$ll}; 
			}
		}
		%lhist = %new_lhist; 
	}
	print STDOUT join("\t", qw/RdLen/, $opts{rd_LenHist_name})."\n"; 
	for my $ll ( $min .. $max ) {
		defined $lhist{$ll} or $lhist{$ll} = 0; 
		print STDOUT join("\t", $ll, $lhist{$ll})."\n"; 
	}
}# rd_LenHist

# Filter reads by read length. 
sub keep_len {
	my ($min_len, $max_len) = &max_min_range( shift ); 

	if ( $opts{paired} ) {
		for (my $i=0; $i<@InFp; $i+=2) {
			my $fh1 = $InFp[$i]; 
			my $fh2 = $InFp[$i+1]; 
			RD: 
			while ( !eof($fh1) and !eof($fh2) ) {
				my $rdRec1 = &get_fq_record($fh1); 
				my $rdRec2 = &get_fq_record($fh2); 
				(my $ss1 = $rdRec1->{seq}) =~ s/\s//g; 
				my $ll1 = length($ss1); 
				$min_len > 0 and $ll1 < $min_len and next RD; 
				$max_len > 0 and $ll1 > $max_len and next RD; 
				(my $ss2 = $rdRec2->{seq}) =~ s/\s//g; 
				my $ll2 = length($ss2); 
				$min_len > 0 and $ll2 < $min_len and next RD; 
				$max_len > 0 and $ll2 > $max_len and next RD; 
				print STDOUT "\@$rdRec1->{id}$rdRec1->{seq}+\n$rdRec1->{qual}\n"; 
				print STDOUT "\@$rdRec2->{id}$rdRec2->{seq}+\n$rdRec2->{qual}\n"; 
			}#End while() RD: 
		}#End for (my $i=0; ...)
	} else {
		for my $fh (@InFp) {
			while ( my $rdRec = &get_fq_record($fh) ) {
				(my $ss = $rdRec->{seq}) =~ s/\s//g; 
				my $ll = length($ss); 
				$min_len > 0 and $ll < $min_len and next; 
				$max_len > 0 and $ll > $max_len and next; 
				print STDOUT "\@$rdRec->{id}$rdRec->{seq}+\n$rdRec->{qual}\n"; 
			}
		}#End for my $fh
	}# if ( $opts{paired} )
}# sub keep_len

sub joinR12 {
	for (my $i=0; $i<@InFp; $i+=2) {
		my $fh1 = $InFp[$i]; 
		my $fh2 = $InFp[$i+1]; 
		RD: 
		while ( !eof($fh1) and !eof($fh2) ) {
			my $rdRec1 = &get_fq_record($fh1); 
			my $rdRec2 = &get_fq_record($fh2); 
			print STDOUT "\@$rdRec1->{id}$rdRec1->{seq}+\n$rdRec1->{qual}\n"; 
			print STDOUT "\@$rdRec2->{id}$rdRec2->{seq}+\n$rdRec2->{qual}\n"; 
		}#End while() RD: 
	}#End for
}#sub joinR12() 



sub sepR12 {
	my $pref = shift; 
	my $of1h = &openFH( "${pref}_R1.fq", '>' ); 
	my $of2h = &openFH( "${pref}_R2.fq", '>' ); 
	for my $fh ( @InFp ) {
		while ( !eof($fh) ) {
			my $rdRec1 = &get_fq_record($fh); 
			my $rdRec2 = &get_fq_record($fh); 
			print {$of1h} "\@$rdRec1->{id}$rdRec1->{seq}+\n$rdRec1->{qual}\n"; 
			print {$of2h} "\@$rdRec2->{id}$rdRec2->{seq}+\n$rdRec2->{qual}\n"; 
		}# while ( !eof() )
	}# for my $fh 
}#sub sepR12() 


# Transvert -phred64 fastq and -phred33 fastq. 
# Another useful method to change quality string is to use 'tr///' command like this (From -phred64 to -phred33): 
# $rdQua =~ tr!>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[ \\ ]^_`abcdefghijklmnopqrstuvwxyz \{ |}~!#####$%&'()*+,-./0123456789:;< \= >?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[ \\ ]^_!;
sub fq2fq {
	my $scale = shift; 
	defined $scale or $scale = 33-64; 
	for my $fh ( @InFp ) {
		while ( my $rdRec = &get_fq_record($fh) ) {
			my $new_qStrR = &char2char( \$rdRec->{qual}, $scale ); 
			print STDOUT "\@$rdRec->{id}$rdRec->{seq}+\n$$new_qStrR\n"; 
		}
	}#
}#sub fq2fq() 


#****************************************************************#
#--------------InnerSubprogram------------Start------------------#
#****************************************************************#

# phred33 to phred64 / phred64 to phred33
# Input : Reference of quality string, and the scale to change char to char. Default is to change phred64 to phred33. 
# Return: Reference of new quality string. 
sub char2char {
	my $raw_qR = shift; 
	my $scale = shift; 
	defined $scale or $scale = 33-64; 
	my @new_ascii; 
	for my $tv ( unpack("C*", $$raw_qR) ) {
		push( @new_ascii, $tv+$scale ); 
	}
	my $back_q = pack("C*", @new_ascii); 
	return(\$back_q); 
}#sub char2char() 


# transform illumina quality character string to phred values. 
# Input : a reference of quality string. And a phred scale used to cut ord() number. 
# Return: an array reference [q1,q2,q3,q4,...]
sub qChar2qValue {
	my $cR = shift; 
	my $scale = shift; 
	defined $scale or $scale = 33; 
	my @back = unpack("C*", $$cR); 
	for my $tv ( @back ) {
		$tv -= $scale; 
	}# for my $tv (@back)
#	pos($$cR) = 0; 
#	while ( $$cR =~ m/(\S)/g ) {
#		push(@back, ord($1)-$scale); 
#	}
	return(\@back); 
}#sub qChar2qValue() 

# Input : an array reference of quality values [q1,q2,q3,...]. And a phred scale used to get character. 
# Return: a reference of character string for the input values. 
sub qValue2qChar {
	my $vR = shift; 
	my $scale = shift; 
	defined $scale or $scale = 33; 
	my $back; 
	for my $tv (@$vR) {
		$back .= chr( $tv + $scale ); 
	}
	return (\$back); 
}

# Read one fastq from file handle once a time. 
# Input : a fastq file handle. 
# Return: A hash reference (\%) with keys: qw/id seq qual/, in which the first '@' in "id" is trimmed, and the quality string is chomped. 
sub get_fq_record {
	my $fh = shift; 
	my %back; 
	eof($fh) and return( undef() ); 
	$back{id} = readline($fh); 
	$back{id} =~ s/^\@// or &stopErr("[Err] The read id line seems wrong:\t$back{id}\n"); 
	$back{'ID'} = $back{'id'}; 
	$opts{rdKey} and $back{id} =~ s/^(\S+)\s.*$/$1/; 
	$back{seq} = readline($fh); 
	readline($fh); 
	$back{qual} = readline($fh); 
	chomp($back{qual});
	return(\%back); 
}

# Guess phred scale from file handle input. 
# Input : File handle of fastq file. 
# Return: Phred scale of fastq quality format. 
sub guessPhredScale {
	defined $opts{phred_scale} and return ($opts{phred_scale}); 
	my $fh = shift; 
	my $not_get = 1; 
	while ($not_get) {
		my $rd = &get_fq_record($fh); 
		chomp($rd->{qual}); 
		pos($rd->{qual}) = 0; 
		while ( $rd->{qual} =~ m/(\S)/g ) {
			my $tQV = ord($1); 
			if ( $tQV <= 58 ) {
				seek( $fh, 0, 0 ); 
				return ( 33 ); # -phred33
			} elsif ( $tQV >= 75 ) {
				seek( $fh, 0, 0 ); 
				return ( 64 ); # -phred64
			}
		}
	}
	&tsmsg("[Err] Failed to guess phred scale for file handle [$fh], use 33 as default.\n"); 
	seek( $fh, 0, 0 ); 
	return( 33 ); # -phred33
}


# input ($seq_ref, $deal_tag); deal_tag : 'r' => reverse, 'c' => complemented, 'rc' => reverse and complemented; Default 'rc';
# no output, edit the input sequence reference.
sub rcSeq {
	my $seq_r = shift;
	my $tag = shift; defined $tag or $tag = 'rc'; # $tag = lc($tag);
	my ($Is_r, $Is_c) = (0)x2;
	$tag =~ /r/i and $Is_r = 1;
	$tag =~ /c/i and $Is_c = 1;
	!$Is_r and !$Is_c and &stopErr( "Wrong Input for function rcSeq! $!\n" );
	$Is_r and $$seq_r = reverse ($$seq_r);
	$Is_c and $$seq_r =~ tr/acgturykmbvdhACGTURYKMBVDHwWsSnN/tgcaayrmkvbhdTGCAAYRMKVBHDwWsSnN/; # edit on 2013-09-11 No difference in result.
	return 0;
}# 2007-9-11 9:46 ÖÆ×÷¶ÔÓ¦·´Ïò»¥²¹ÐòÁÐ;
#        a       a; adenine
#        c       c; cytosine
#        g       g; guanine
#        t       t; thymine in DNA; uracil in RNA
#        m       a or c
#        k       g or t
#        r       a or g
#        y       c or t
#        w       a or t
#        s       c or g
#        v       a or c or g; not t
#        b       c or g or t; not a
#        h       a or c or t; not g
#        d       a or g or t; not c
#        n       a or c or g or t

# Input : "min_len-max_len"
# Return: (min_len, max_len)
sub max_min_range {
	my $kplen = shift; 
	my ($min_len, $max_len) = (-1, -1); 
	if ( $kplen =~ m!^\s*(\d+)\-(\d+)\s*$! ) {
		($min_len, $max_len) = ($1, $2);
	} elsif ( $kplen =~ m!^\s*\-(\d+)\s*$! ) {
		$max_len = $1;
	} elsif ( $kplen =~ m!^\s*(\d+)\-\s*$! ) {
		$min_len = $1;
	} else {
		&stopErr( "[Err] Failed to parse option [-keep_len $kplen]\n" );
	}
	return($min_len, $max_len); 
}# sub max_min_range

# Input : string to give a set of ranges. 
# Return: A reference of array whose element is array of [start, end]
sub frag_regions {
	my $frag_str = shift; 
	( defined $frag_str and $frag_str ne '') or &stopErr("[Err] No value given to \$frag_str.\n"); 
	my @back; 
	for ( split(/:/, $frag_str) ) {
		s/^\s+//; s/\s+$//; 
		if ( m/^(\-?(?:\d+)?)-(\-?(?:\d+)?)$/ ) {
			my ($ss, $ee) = ($1,$2); 
			push( @back, [ ($ss) ? $ss : 1, ($ee) ? $ee : 'end' ] ); 
			(defined $ss and defined $ee and $ss > $ee and $ee > 0) and &stopErr("[Err] Bad fragment region [$frag_str].\n"); 
		} else {
			&stopErr( "[Err] fragment region [$frag_str] not accepted.\n" ); 
		}
	}
	return(\@back); 
}#sub frag_regions() 


