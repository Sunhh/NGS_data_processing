package SeqAlnSunhh; 

#Author: Honghe Sun, hs738@cornell.edu, 2014-07-15
# 2014-07-15 In function olap_e2e_A2B(), "N" is always treated as "mismatch"; 

use strict; 
use warnings; 
use LogInforSunhh; 
use Exporter qw(import); 
use mathSunhh; # For usage of mathSunhh->_setHashFromArr(); 
my $ms=mathSunhh->new(); 

our @EXPORT = qw(olap_e2e_A2B); 
our @EXPORT_OK; 

############################################################
#  Basic settings. 
############################################################

## For olap_e2e_A2B() sub-routine. 
# define the scores for match, mismatch and indel:
my %para; 
$para{match_score} = 2; 
$para{mismatch_score} = -1; 
$para{gapopen_score} = 0; 
$para{gapextend_score} = -2; 
$para{opt_tail} = 0; 
$para{opt_head} = 0; 



############################################################
#  Methods
############################################################
sub new {
	my $class = shift; 
	my $self = {}; 
	bless $self, $class; 

	$self->_initialize(@_); 

	return $self; 
}
sub _initialize {
	my $self = shift; 
	my %parm = $ms->_setHashFromArr(@_); 
	return; 
}

=head2 aln_tophat2( inFq1=>[$PE_Fq1, $PE_Fq2, ...], inFq2=>[$PE_Fq2, $PE_Fq2, ...], 
 index=>$bowtie_index, outDir=>'./tophat_out', para_tophat=>'', printCmd=>0 )

Description   : 

Input         : 
 inFq1      : Requried. For SE fastq, all files should be given to inFq1. 
 inFq2      : Provided along with inFq1 for PE fastq. 
 index      : bowtie2 index built with bowtie2-build
 para_tophat: For example: '-p $cpuN --library-type=fr-firststrand --read-mismatches 1 --splice-mismatches 0 --min-intron-length 30'
 outDir     : The directory recording all output. We will use '-o outDir' in '-para_tophat' if it is given. 
 printCmd   : Only print commands if given. 
=cut
sub aln_tophat2 {
	my $self = shift; 
	my %parm = $ms->_setHashFromArr(@_); 
	
	# Find tophat2 
	$parm{'exe_tophat'} //= 'tophat'; 
	# Check tophat version: 
	my $tophat_version = `$parm{'exe_tophat'} --version`; 
	chomp($tophat_version); 
	$tophat_version =~ m/^\s*tophat\s*v2\./i or &tsmsg("[Err]The tophat version is [$tophat_version] instead of v2, which may cause problem!\n"); 
	
	# Set Other parameters. 
	$parm{'para_tophat'} //= ''; 
	$parm{'printCmd'} //= 0; 
	$parm{'outDir'} //= './tophat_out'; 
	if ( $parm{'para_tophat'} =~ s!(?:\-o|\-\-output\-dir)\s+(\S+)(?:\s*$|\s+\-)!! ) {
		my $replace_dir = $1; 
		&tsmsg("[Wrn]Replace outDir from [$parm{'outDir'}] to [$replace_dir]\n"); 
		$parm{'outDir'} = $replace_dir; 
	}
	
	# Separate and group PE/SE fastq files. 
	$parm{'inFq1'} //= []; 
	$parm{'inFq2'} //= []; 
	my $long_idx = $ms->max( $#{$parm{'inFq1'}}, $#{$parm{'inFq2'}} ); 
	my (@inFqPE1, @inFqPE2, @inFqSE); 
	for (my $i=0; $i<=$long_idx; $i++) {
		if ( defined $parm{'inFq1'}[$i] and $parm{'inFq1'}[$i] ne '' ) {
			if ( defined $parm{'inFq2'}[$i] and $parm{'inFq2'}[$i] ne '' ) {
				push(@inFqPE1, $parm{'inFq1'}[$i]); 
				push(@inFqPE2, $parm{'inFq2'}[$i]); 
			} else {
				push(@inFqSE, $parm{'inFq1'}[$i]); 
			}
		} elsif ( defined $parm{'inFq2'}[$i] and $parm{'inFq2'}[$i] ne '' ) {
			push(@inFqSE, $parm{'inFq2'}[$i]); 
		} else {
			# Faint. 
		}
	}
	
	# Check bowtie2-index 
	$parm{'printCmd'} or $self->chk_index( $parm{'index'}, 'type'=>'bowtie2' ); 
	
	# Run tophat2 
	my ($pe_str1, $pe_str2); 
	$pe_str1 = join(",", @inFqPE1, @inFqSE); 
	$pe_str2 = join(",", @inFqPE2); 
	$pe_str1 //= ''; 
	$pe_str2 //= ''; 
	&exeCmd_1cmd("$parm{'exe_tophat'} $parm{'para_tophat'} --output-dir $parm{'outDir'} $parm{'index'} $pe_str1 $pe_str2", $parm{'printCmd'}); 
	
	return; 
}# aln_tophat2

=head2 chk_index( $index_prefix, 'type'=>'bowtie2' )
Description   : Check database index files with given 'index_prefix' and database 'type'
=cut
sub chk_index {
	my $self = shift; 
	my $pref = shift; 
	my %parm = $ms->_setHashFromArr(@_); 
	
	$parm{'type'} //= 'bowtie2'; 
	$parm{'type'} = lc($parm{'type'}); 
	if ( $parm{'type'} eq 'bowtie2' ) {
		for my $suff (qw/.1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2/) {
			-e "${pref}${suff}" or do { &tsmsg("[Err] index file [${pref}${suff}] not found.\n"); return 1; }; 
		}
		return 0; 
	} else {
		&tsmsg("[Err]Unknown database type [$parm{'type'}] for index [$pref]\n"); 
		return 1; 
	}
	# Why here!
	return 1; 
}# sub chk_index() 


=head2 bwaAln( inFq1=>$inFq1, aln_type=>'PE', %parameters_for_bwaPE_or_bwaSE_function )

Description   : Invoke bwaPE()/bwaSE() to run bwa alignment. 

Input         : 
  aln_type : Default 'PE' = paired alignment. Can be 'SE' = single alignments. 

=cut

sub bwaAln {
	my $self = shift; 
	my %parm = $ms->_setHashFromArr(@_); 
	$parm{'aln_type'} //= 'PE'; 
	if ( $parm{'aln_type'} eq 'PE' ) {
		$self->bwaPE(%parm); 
	} elsif ( $parm{'aln_type'} eq 'SE' ) {
		$self->bwaSE(%parm); 
	} else {
		&tsmsg("[Err] Unknown aln_type [$parm{'aln_type'}]\n"); 
	}
}# bwaAln() 

=head2 bwaPE( inFq1=>$inFq1, inFq2=>$inFq2, step2do=>['all'], db=>$bwa_aln_db, dbTag=>'toDB', 
  oBamPre=>'bwaPEOutPre', 
  exe_bwa=>'bwa', exe_samtools=>'samtools', 
  para_aln=>'', para_sampe=>'', 
  para_bamSort=>'', 
  para_sam2bam=>'', 
  printCmd=>0, 
)

Description   : Please note that I don't check the file existence! 

Input         : 
  step2do : Could be [ qw/all aln_r1 aln_r2 sampe bam_sort bam_index rm_sai rm_rawbam sampe2bam sampe2sam sam2bam rm_rawsam/ ]
            Here 'all' means qw/aln_r1 aln_r2 sampe bam_sort bam_index rm_sai rm_rawbam/
  printCmd: 0 (default) indicates we carry out each command. 1 indicates we only print command with tag [CMD_print]. 

=cut
sub bwaPE {
	my $self = shift; 
	my %parm = $ms->_setHashFromArr(@_); 
	( ( defined $parm{'inFq1'} or defined $parm{'inFq2'} ) and defined $parm{'db'} ) or &stopErr("[Err] Three paramters must be defined: inFq1/inFq2/db in function bwaPE()\n"); 
	$parm{'step2do'} //= ['all']; 
	$parm{'dbTag'} //= 'toDB'; 
	$parm{'exe_bwa'} //= 'bwa'; 
	$parm{'exe_samtools'} //= 'samtools'; 
	$parm{'para_aln'} //= ''; # Recommanded : -t $cpuN -n 0.03 -o 1 -e 2
	$parm{'para_sampe'} //= ''; # Recommanded : 
	$parm{'para_bamSort'} //= ''; # Recommanded : -@ $cpuN -m $mem_limit
	$parm{'para_sam2bam'} //= ''; 
	$parm{'oBamPre'} //= 'bwaPEOutPre'; 
	$parm{'printCmd'} //= 0; 

	# Setting file names. 
	my $inFq1 = $parm{'inFq1'}; 
	my $inFq2 = $parm{'inFq2'}; 
	my $oBamPre = $parm{'oBamPre'}; 
	my $db = $parm{'db'}; 
	my $dbTag = $parm{'dbTag'}; 
	my $oSai1 = "${inFq1}.${dbTag}.sai"; 
	my $oSai2 = "${inFq2}.${dbTag}.sai"; 

	for my $t_step ( @{$parm{'step2do'}} ) {
		my @steps = (lc($t_step)); 
		$t_step eq 'all' and @steps = (qw/aln_r1 aln_r2 sampe bam_sort bam_index rm_sai rm_rawbam/); 
		for my $c_step ( @steps ) {
			if      ( $c_step eq 'aln_r1' ) {
				&exeCmd_1cmd("$parm{'exe_bwa'} aln $parm{'para_aln'} -f $oSai1 $db $inFq1", $parm{'printCmd'}); 
			} elsif ( $c_step eq 'aln_r2' ) {
				&exeCmd_1cmd("$parm{'exe_bwa'} aln $parm{'para_aln'} -f $oSai2 $db $inFq2", $parm{'printCmd'}); 
			} elsif ( $c_step eq 'sampe' ) { 
				&exeCmd_1cmd("$parm{'exe_bwa'} sampe $parm{'para_sampe'} $db $oSai1 $oSai2 $inFq1 $inFq2 | $parm{'exe_samtools'} view $parm{'para_sam2bam'} -bSh -o $oBamPre.bam -", $parm{'printCmd'}); 
			} elsif ( $c_step eq 'sampe2bam' ) {
				&exeCmd_1cmd("$parm{'exe_bwa'} sampe $parm{'para_sampe'} $db $oSai1 $oSai2 $inFq1 $inFq2 | $parm{'exe_samtools'} view $parm{'para_sam2bam'} -bSh -o $oBamPre.bam -", $parm{'printCmd'}); 
			} elsif ( $c_step eq 'sampe2sam' ) {
				&exeCmd_1cmd("$parm{'exe_bwa'} sampe $parm{'para_sampe'} $db $oSai1 $oSai2 $inFq1 $inFq2 > $oBamPre.sam", $parm{'printCmd'}); 
			} elsif ( $c_step eq 'sam2bam' ) {
				&exeCmd_1cmd("$parm{'exe_samtools'} view $parm{'para_sam2bam'} -bSh -o $oBamPre.bam $oBamPre.sam", $parm{'printCmd'}); 
			} elsif ( $c_step eq 'bam_sort' ) {
				&exeCmd_1cmd("$parm{'exe_samtools'} sort $parm{'para_bamSort'} $oBamPre.bam $oBamPre.srt", $parm{'printCmd'}); 
			} elsif ( $c_step eq 'bam_index' ) {
				&exeCmd_1cmd("$parm{'exe_samtools'} index $oBamPre.srt.bam", $parm{'printCmd'}); 
			} elsif ( $c_step eq 'rm_sai' ) {
				my @file_list; 
				if ( $parm{'printCmd'} ) {
					@file_list = ( $oSai1, $oSai2 ); 
				} else {
					@file_list = grep { -f $_ } ( $oSai1, $oSai2 ); 
				}
				$#file_list > -1 and &exeCmd_1cmd("rm " . join(" ", @file_list), $parm{'printCmd'}); 
			} elsif ( $c_step eq 'rm_rawbam' ) {
				if ( $parm{'printCmd'} ) {
					&exeCmd_1cmd("rm $oBamPre.bam", $parm{'printCmd'}); 
				} else {
					-f "$oBamPre.bam" and &exeCmd_1cmd("rm $oBamPre.bam", $parm{'printCmd'}); 
				}
			} elsif ( $c_step eq 'rm_rawsam' ) {
				if ( $parm{'printCmd'} ) {
					&exeCmd_1cmd("rm $oBamPre.sam", $parm{'printCmd'}); 
				} else {
					-f "$oBamPre.sam" and &exeCmd_1cmd("rm $oBamPre.sam", $parm{'printCmd'}); 
				}
			} else {
				&tsmsg("[Err] Unknown step [$c_step]\n"); 
			}
		}
	}

	return; 
}# End bwaPE() 

=head2 bwaSE( inFq1=>$inFq1, , step2do=>['all'], db=>$bwa_aln_db, dbTag=>'toDB', 
  oBamPre=>'bwaSEOutPre', 
  exe_bwa=>'bwa', exe_samtools=>'samtools', 
  para_aln=>'', para_samse=>'', 
  para_bamSort=>'', 
  para_sam2bam=>'', 
  printCmd=>0, 
)
Description   : Please note that I don't check the file existence! 

Input         : 
  step2do : Could be [ qw/all aln_r1 samse bam_sort bam_index rm_sai rm_rawbam samse2sam sam2bam rm_rawsam/ ]
            Here 'all' means qw/aln_r1 samse bam_sort bam_index rm_sai rm_rawbam/
  printCmd: 0 (default) indicates we carry out each command. 1 indicates we only print command with tag [CMD_print]. 

=cut
sub bwaSE {
	my $self = shift; 
	my %parm = $ms->_setHashFromArr(@_); 
	( ( defined $parm{'inFq1'} ) and defined $parm{'db'} ) or &stopErr("[Err] Three paramters must be defined: inFq1/db in function bwaSE()\n"); 
	$parm{'step2do'} //= ['all']; 
	$parm{'dbTag'} //= 'toDB'; 
	$parm{'exe_bwa'} //= 'bwa'; 
	$parm{'exe_samtools'} //= 'samtools'; 
	$parm{'para_aln'} //= ''; 
	$parm{'para_samse'} //= ''; 
	$parm{'para_bamSort'} //= ''; 
	$parm{'para_sam2bam'} //= ''; 
	$parm{'oBamPre'} //= 'bwaSEOutPre'; 
	$parm{'printCmd'} //= 0; 

	# Setting file names. 
	my $inFq1 = $parm{'inFq1'}; 
	my $oBamPre = $parm{'oBamPre'}; 
	my $db = $parm{'db'}; 
	my $dbTag = $parm{'dbTag'}; 
	my $oSai1 = "${inFq1}.${dbTag}.sai"; 

	for my $t_step ( @{$parm{'step2do'}} ) {
		my @steps = (lc($t_step)); 
		$t_step eq 'all' and @steps = (qw/aln_r1 samse bam_sort bam_index rm_sai rm_rawbam/); 
		for my $c_step ( @steps ) {
			if      ( $c_step eq 'aln_r1' ) {
				&exeCmd_1cmd("$parm{'exe_bwa'} aln $parm{'para_aln'} -f $oSai1 $db $inFq1", $parm{'printCmd'}); 
			} elsif ( $c_step eq 'samse' ) { 
				&exeCmd_1cmd("$parm{'exe_bwa'} samse $parm{'para_samse'} $db $oSai1 $inFq1 | $parm{'exe_samtools'} view $parm{'para_sam2bam'} -bSh -o $oBamPre.bam -", $parm{'printCmd'}); 
			} elsif ( $c_step eq 'samse2bam' ) {
				&exeCmd_1cmd("$parm{'exe_bwa'} samse $parm{'para_samse'} $db $oSai1 $inFq1 | $parm{'exe_samtools'} view $parm{'para_sam2bam'} -bSh -o $oBamPre.bam -", $parm{'printCmd'}); 
			} elsif ( $c_step eq 'samse2sam' ) {
				&exeCmd_1cmd("$parm{'exe_bwa'} samse $parm{'para_samse'} $db $oSai1 $inFq1 > $oBamPre.sam", $parm{'printCmd'}); 
			} elsif ( $c_step eq 'sam2bam' ) {
				&exeCmd_1cmd("$parm{'exe_samtools'} view $parm{'para_sam2bam'} -bSh -o $oBamPre.bam $oBamPre.sam", $parm{'printCmd'}); 
			} elsif ( $c_step eq 'bam_sort' ) {
				&exeCmd_1cmd("$parm{'exe_samtools'} sort $parm{'para_bamSort'} $oBamPre.bam $oBamPre.srt", $parm{'printCmd'}); 
			} elsif ( $c_step eq 'bam_index' ) {
				&exeCmd_1cmd("$parm{'exe_samtools'} index $oBamPre.srt.bam", $parm{'printCmd'}); 
			} elsif ( $c_step eq 'rm_sai' ) {
				my @file_list; 
				if ( $parm{'printCmd'} ) {
					@file_list = ( $oSai1 ); 
				} else {
					@file_list = grep { -f $_ } ( $oSai1 ); 
				}
				$#file_list > -1 and &exeCmd_1cmd("rm " . join(" ", @file_list), $parm{'printCmd'}); 
			} elsif ( $c_step eq 'rm_rawbam' ) {
				if ( $parm{'printCmd'} ) {
					&exeCmd_1cmd("rm $oBamPre.bam", $parm{'printCmd'}); 
				} else {
					-f "$oBamPre.bam" and &exeCmd_1cmd("rm $oBamPre.bam", $parm{'printCmd'}); 
				}
			} elsif ( $c_step eq 'rm_rawsam' ) {
				if ( $parm{'printCmd'} ) {
					&exeCmd_1cmd("rm $oBamPre.sam", $parm{'printCmd'}); 
				} else {
					-f "$oBamPre.sam" and &exeCmd_1cmd("rm $oBamPre.sam", $parm{'printCmd'}); 
				}
			} else {
				&tsmsg("[Err] Unknown step [$c_step]\n"); 
			}
		}
	}

}# End bwaSE() 



############################################################
#  Sub-routines. 
############################################################


=head1 olap_e2e_A2B( $sequenceA, $sequenceB, [, {%para}] )

Description : This function is used to align two sequences from end to end. 
              Do not distinguish "N" or other special words. 

Output      : { qw/(start|end)(A|B) aln_len aln_ident count seqA_aln seqB_aln seqC_aln/ => values }
                startA, endA : start and end positions of seqA aligned; 
                startB, endB : start and end positions of seqB aligned; 
                aln_len      : length of aligned region (including gap)
                aln_ident    : (length of match bases) / aln_len
                seqA_aln     : alignment in seqA
                seqB_aln     : alignemtn in seqB
                seqC_aln     : Consensus signs between seqA_aln and seqB_aln. "." - same, "D" - diff, "-" - gap. 

=cut
sub olap_e2e_A2B {
	my $seqA = shift; 
	unless ( ref($seqA) eq 'SCALAR' ) {
		ref($seqA) eq 'SeqAlnSunhh' or &stopErr("[Err] olap_e2e_A2B 1st-input should be a scalar string standing for sequence A.\n"); 
		$seqA = shift; 
	}
	my ($seqB, $ref_para) = @_; 
	my %loc_para = %para; 
	if (defined $ref_para and ref($ref_para) eq 'HASH') {
		for my $tk ( keys %$ref_para ) {
			$loc_para{$tk} = $ref_para->{$tk}; 
		}
	}
	
	# return back values; 
	my ($psA, $peA, $psB, $peB); # (startA, endA, startB, endB) of overlap region. 
	my $seqA_aln = ''; # Aligned sequence from seqA 
	my $seqB_aln = ''; # Aligned sequence from seqB 
	my ($aln_len, $aln_ident) = (0, 0); 

	# Caculating. 
	$seqA = uc($seqA); 
	$seqB = uc($seqB); 
	my $lenA = length($seqA); 
	my $lenB = length($seqB); 
	my @baseA = ( $seqA =~ m/(.)/g ); 
	my @baseB = ( $seqB =~ m/(.)/g ); 
	my @smat;               # array matrix storing the dynamic programming score matrix; 
	my @traceback;          # array storing the traceback directions. 
	
	# Define two-dimensional matrices "@smat" and "@traceback" with ($lenA + 1) rows, and ($lenB + 1) columns, and initialize it with 0s. 
	for (my $i=0; $i<=$lenA; $i++) {
		for (my $j=0; $j<=$lenB; $j++) {
			$smat[$i][$j]      = 0; 
			$traceback[$i][$j] = 0; 
		}
	}
	
	# Setup boundary of traceback array, used to stop searching. 
	# Put "."s in the first row of matrix @traceback; 
	for (my $j=0; $j<=$lenB; $j++) {
		$traceback[0][$j] = "."; 
	}
	# Put "."s in the first column of matrix @traceback; 
	for (my $i=0; $i<=$lenA; $i++) {
		$traceback[$i][0] = "."; 
	}
	
	# Dynamic programming recursion
	for (my $i=1; $i<=$lenA; $i++) {
		# base ($baseA[$i-1]) position $i in $seqA (seqA goes down the first column) 
		for (my $j=1; $j<=$lenB; $j++) {
			# base ($baseB[$j-1]) position $j in $seqB (seqB goes across the first row)
			my ($score , $diag, $up, $left, $max) ; # scores used to construct matrix; 
			
			# Find the value to put into $smat[$i][$j]
			# (1) the first possibility is to take the diagonal element plus a match/mismatch score; 
			if ( $baseA[$i-1] eq $baseB[$j-1] and $baseA[$i-1] ne "N" ) { $score = $loc_para{match_score};    } 
			else                                { $score = $loc_para{mismatch_score}; }
			$diag = $smat[$i-1][$j-1]+$score; 
			# (2) the second possibility is to take the element above plus a gap score; 
			if ( $traceback[$i-1][$j  ] eq '|' ) { $up = $smat[$i-1][$j  ] + $loc_para{gapextend_score};                        } 
			else                                 { $up = $smat[$i-1][$j  ] + $loc_para{gapextend_score} + $loc_para{gapopen_score}; }
			# (3) the third possibility is to take the element on the left plus a gap score; 
			if ( $traceback[$i  ][$j-1] eq '-' ) { $left = $smat[$i  ][$j-1] + $loc_para{gapextend_score};                        }
			else                                 { $left = $smat[$i  ][$j-1] + $loc_para{gapextend_score} + $loc_para{gapopen_score}; }
			# record which of the three possibilities was highest, into @traceback: 
			if    ( $diag >  $up   && $diag >  $left ) { $traceback[$i][$j] = '>'; $max = $diag; }
			elsif ( $up   >  $diag && $up   >  $left ) { $traceback[$i][$j] = '|'; $max = $up;   }
			elsif ( $left >  $diag && $left >  $up   ) { $traceback[$i][$j] = '-'; $max = $left; }
			elsif ( $left == $diag && $up   == $diag ) { $traceback[$i][$j] = '*'; $max = $diag; }
			elsif ( $up   == $left && $up   >  $diag && $left >  $diag ) { $traceback[$i][$j] = 'L'; $max = $up;   }
			elsif ( $up   == $diag && $up   >  $left && $diag >  $left ) { $traceback[$i][$j] = 'V'; $max = $up;   }
			elsif ( $diag == $left && $diag >  $up   && $left >  $up   ) { $traceback[$i][$j] = 'Z'; $max = $diag; }
			else { &stopErr( "Why we come here!\n" );  }
			# record the highest score in @smat; 
			$smat[$i][$j] = $max; 
		}
	}
	
	# The best overlap is the given by the largest number in the bottom-most row of the matrix @smat; 
	# Because we are doing A-2-B overlap alignment; 
	my $best_x = -1; # Position aligned in seqA 
	my $best_y = -1; # Position aligned in seqB 
	my $best_score = -1e15; 
	
	# Here we only check the end of seqA, and we want to get the left-most alignment for A-2-B end-to-end overlap 
	# So we check only the bottom-most row from left to right. 
	if ($loc_para{opt_tail} > 0) {
		for (my $i=$lenA; $i>=0 && $i >= $lenA-$loc_para{opt_tail} ; $i--) {
			for (my $j=0; $j<=$lenB; $j++) {
				if ( $smat[$i][$j] > $best_score ) { $best_x = $i; $best_y = $j; $best_score = $smat[$i][$j]; } 
			}
		}
	}else{
		for (my $j=0; $j<=$lenB; $j++) { 
			# If we want right-most alignment on seqB, we should use "for (my $j=$lenB; $j>=0; $j--)"
			if ( $smat[$lenA][$j] > $best_score ) { $best_x = $lenA; $best_y = $j; $best_score = $smat[$lenA][$j]; } 
		}
	}
	
	# Record the boundary of best overlap region. 
	$peA = $best_x; 
	$peB = $best_y; 
	
	## record the traceback as '#': 
	# construct the best alignment. 
	my $tracebackvalue = $traceback[$best_x][$best_y]; 
	my $prev_score = $smat[$best_x][$best_y]; 
	my %count = qw(same 0 diff 0 gapA 0 gapB 0); 
	my $seqC_aln = ''; 
	while ( $tracebackvalue ne '.' ) {
		$traceback[$best_x][$best_y] = "#$traceback[$best_x][$best_y]"; 
		my ($prev_x, $prev_y) = ($best_x, $best_y); 
		my ($seqA_letter, $seqB_letter, $seqC_letter, $cKey); 
		if ( $tracebackvalue eq '>' or $tracebackvalue eq '*' or $tracebackvalue eq 'V' or $tracebackvalue eq 'Z' ) {
			$seqA_letter = $baseA[$best_x-1]; 
			$seqB_letter = $baseB[$best_y-1]; 
			if ( $seqA_letter eq $seqB_letter ) {
				$cKey = 'same'; 
				$seqC_letter = '.'; 
			} else {
				$cKey = 'diff'; 
				$seqC_letter = 'D'; 
			}
			$best_x --; 
			$best_y --; 
		} elsif ( $tracebackvalue eq '-' or $tracebackvalue eq 'L' ) {
			$seqA_letter = '-'; 
			$seqB_letter = $baseB[$best_y-1]; 
			$cKey = 'gapA'; 
			$seqC_letter = '-'; 
			$best_y --; 
		} elsif ( $tracebackvalue eq '|' ) {
			$seqA_letter = $baseA[$best_x-1]; 
			$seqB_letter = '-'; 
			$cKey = 'gapB'; 
			$seqC_letter = '-'; 
			$best_x --; 
		} else {
			&tsmsg( "[Err] Unknown tracebackvalue [$tracebackvalue]\n" ); 
			exit 1; 
		}
		$seqA_aln = $seqA_letter . $seqA_aln; 
		$seqB_aln = $seqB_letter . $seqB_aln; 
		$seqC_aln = $seqC_letter . $seqC_aln; 
		$count{$cKey} ++; 
		$aln_len ++; 
		if ( $loc_para{opt_head} > 0 ) {
			if ( $best_y <= $loc_para{opt_head} and $prev_score < $smat[$best_x][$best_y] ) {
				# ($best_x, $best_y) = ($prev_x, $prev_y); 
				last; 
			}
		}
		# Find the new traceback value; 
		$tracebackvalue = $traceback[$best_x][$best_y]; 
		$prev_score = $smat[$best_x][$best_y]; 
	}
	$psA = $best_x+1; 
	$psB = $best_y+1; 
	
	## Prepare to output score matrix (@smat) and traceback matrix (@traceback)
	## Put $seqB in the first row of matrix @smat; 
	for (my $j=1; $j<=$lenB; $j++) {
		$smat[0][$j] = $baseB[$j-1]; 
	}
	## Put $seqA in the first column of matrix @smat; 
	for (my $i=1; $i<=$lenA; $i++) {
		$smat[$i][0] = $baseA[$i-1]; 
	}
#	## Print out the result. 
#	print "Overlap matrix, with traceback shown as #:\n"; 
#	my $smatvalue; 
#	$tracebackvalue = undef(); 
#	for (my $i=0; $i<=$lenA; $i++) {
#		for (my $j=0; $j<=$lenB; $j++) {
#			$smatvalue = $smat[$i][$j]; 
#			$tracebackvalue = $traceback[$i][$j]; 
#			print "$smatvalue $tracebackvalue\t"; 
#		}
#		print "\n"; 
#	}
#	print "\n"; 
	
#	# Print out the overlap alignment; 
#	print "Overlap alignment (score = $best_score):\n"; 
#	print "A: $seqA_aln\n"; 
#	print "B: $seqB_aln\n"; 
#	print "C: $seqC_aln\n"; 
	
	$aln_ident = ( $aln_len == 0 ) ? -1 : ($count{same})/$aln_len ; 
	my %back; 
	$back{startA} = $psA; $back{endA} = $peA; 
	$back{startB} = $psB; $back{endB} = $peB; 
	$back{aln_len} = $aln_len; $back{aln_ident} = $aln_ident; 
	$back{aln_score} = $best_score; 
	$back{count} = \%count; 
	$back{seqA_aln} = $seqA_aln; $back{seqB_aln} = $seqB_aln; 
	$back{seqC_aln} = $seqC_aln; 
	return \%back; 
}#End sub olap_e2e_A2B

1; # It is important to include this line. 

