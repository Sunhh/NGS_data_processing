package SeqAlnSunhh; 

#Author: Honghe Sun, hs738@cornell.edu, 2014-07-15
# 2014-07-15 In function olap_e2e_A2B(), "N" is always treated as "mismatch"; 
# 2015-07-01 Add sam_flag_infor() function. 

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

my %infor_flag; 
$infor_flag{ 0 } = [qw/0     0x0001  p/, "the read is paired in sequencing"];
$infor_flag{ 1 } = [qw/1     0x0002  P/, "the read is mapped in a proper pair"];
$infor_flag{ 2 } = [qw/2     0x0004  u/, "the query sequence itself is unmapped"];
$infor_flag{ 3 } = [qw/3     0x0008  U/, "the mate is unmapped"];
$infor_flag{ 4 } = [qw/4     0x0010  r/, "strand of the query (1 for reverse)"];
$infor_flag{ 5 } = [qw/5     0x0020  R/, "strand of the mate (1 for reverse)"];
$infor_flag{ 6 } = [qw/6     0x0040  1/, "the read is the first read in a pair"];
$infor_flag{ 7 } = [qw/7     0x0080  2/, "the read is the second read in a pair"];
$infor_flag{ 8 } = [qw/8     0x0100  s/, "the alignment is not primary"];
$infor_flag{ 9 } = [qw/9     0x0200  f/, "the read fails platform/vendor quality checks"];
$infor_flag{ 10} = [qw/10    0x0400  d/, "the read is either a PCR or an optical duplicate"];
$infor_flag{ 11} = [qw/11    0x0800  d/, "the alignment is supplementary alignment"]; 

my %flag_grp; 

# _chk_sam_flag( $flag_key, $flag_to_chk )
# This is used to check if flag_to_chk_value is defined in group 'flag_key' ; 
# Returns 1-for yes or 0-for no. 
sub _chk_sam_flag {
	my ($flag_key, $flag_to_chk, $debug) = @_; 
	$debug //= 1; 
	$debug and !(defined $flag_grp{$flag_key} ) and &tsmsg("[Debug] Generating sam_flag for [$flag_key]\n"); 
	if ($flag_key eq 'hDiff_Forward') {
		$flag_grp{$flag_key} //= &mk_flag( 'keep'=>'0=1,2=0,3=0,4=0,5=1' , 'drop'=>'' ); 
	} elsif ($flag_key eq 'hDiff_Reverse') {
		$flag_grp{$flag_key} //= &mk_flag( 'keep'=>'0=1,2=0,3=0,4=1,5=0' , 'drop'=>'' ); 
	} elsif ($flag_key eq 'hDiff_Pair') {
		$flag_grp{$flag_key} //= &mk_flag( 'keep'=>'0=1,2=0,3=0,4=0,5=1;0=1,2=0,3=0,4=1,5=0', 'drop'=>'' ); 
	} elsif ($flag_key eq 'pair_aligned') {
		$flag_grp{$flag_key} //= &mk_flag( 'keep'=>'0=1,2=0,3=0', 'drop'=>'' ); 
	} elsif ($flag_key eq 'read_aligned') {
		$flag_grp{$flag_key} //= &mk_flag( 'keep'=>'2=0', 'drop'=>'' ); 
	} elsif ($flag_key eq 'forward_aligned') {
		$flag_grp{$flag_key} //= &mk_flag( 'keep'=>'2=0,4=0', 'drop'=>'' ); 
	} elsif ($flag_key eq 'reverse_aligned') {
		$flag_grp{$flag_key} //= &mk_flag( 'keep'=>'2=0,4=1', 'drop'=>'' ); 
	} elsif ($flag_key =~ m/^base_(\d+)$/) {
		( $1 >= 0 and $1 <= 10 ) or &stopErr("[Err] flag_key [$flag_key] not accepted.\n"); 
		$flag_grp{$flag_key} //= &mk_flag( 'keep'=>"$1=1", 'drop'=>'' ); 
	} else {
		&stopErr("[Err] flag_key [$flag_key] not accepted.\n"); 
	} 
	my $back = ( defined $flag_grp{$flag_key}{$flag_to_chk} ) ? 1 : 0 ; 
	return $back; 
}

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

=head1 not_uniqBest( \@sam_line_array )

Function : Judge if a SAM alignment is a not-good alignment: (any of the following applies)
            Rule 1 : Contain 'XT:A:R' ; 
            Rule 2 : Have 'XA:Z:*' tag and NM in it is not larger than raw NM:i_value;

Return   : (1/0)

=cut
sub not_uniqBest {
	my $ar = shift; # [@sam_line]
	my ($nm, $xt, $xa); 
	for (my $i=11; $i<@$ar; $i++) {
		local $_ = $ar->[$i];
		if (m/^XT:A:(\S+)$/) {
			defined $xt and die "repeat XT:A in : @$ar\n";
			$xt = $1;
		} elsif (m/^NM:i:(\d+)$/) {
			defined $nm and die "repeat NM:i in : @$ar\n";
			$nm = $1;
		} elsif (m/^XA:Z:(\S+)$/) {
			defined $xa and die "repeat XA:Z in : @$ar\n";
			my $str = $1;
			# XA:Z:scaffold116_cov122,+373278,113M,0;scaffold200_cov109,-2398803,113M,1;
			for my $ts0 (split(/;/, $str)) {
				$ts0 =~ m/^\s*$/ and next;
				$ts0 =~ m/^(\S+),([+-]?\d+),([\d\w]+),(\d+)$/ or die "Failed for XA:Z : [$ts0]\n";
				my ($chr, $pos, $cigar, $nm_1) = ($1, $2, $3, $4);
				$xa //= $nm_1;
				$xa > $nm_1 and $xa = $nm_1;
			}
		} else {
			;
		}
	}
	$nm //= 0;
	$xt //= 'U';
	$xa //= 99999;
	$xa <= $nm and return 1;
	$xt eq 'R' and return 1;
	return 0;
}# sub not_uniqBest() 



=head1 mk_flag( 'keep'=>'', 'drop'=>'' )

Function: 

Return    : \%good_flag_list

%good_flag_list == ( $flag_to_keep1 => 1 , $flag_to_keep2 => 1, ... )

Example: 
 -anyEnd_pair    'keep'=>'0=1' , 'drop'=>'2=1,3=1' ; 
 -anyEnd_info    'keep'=>''    , 'drop'=>'2=1,3=1' ; 
 -bothEnd_pair   'keep'=>'0=1,2=0,3=0' , 'drop'=>'' ; 
 -h2diff_pair    'keep'=>'0=1,2=0,3=0,4=0,5=1;0=1,2=0,3=0,4=1,5=0' , 'drop'=>'' ; 
 -h2diff_F       'keep'=>'0=1,2=0,3=0,4=0,5=1' , 'drop'=>'' ; 
 -onlyMapRd      'keep'=>'', 'drop'=>'2=1' ; 
 -plusRd         'keep'=>'' , 'drop'=>'2=0,4=0'; 
 -minusRd        'keep'=>'2=0,4=1' , 'drop'=>''; 

=cut
sub mk_flag {
	my %parm = @_; 
	$parm{'keep'} //= ''; 
	$parm{'drop'} //= ''; 
	
	my %flag = %{ &_sam_flag_href() }; 
	my ($have_drop, $have_keep) = (0,0); 
	$parm{'keep'} eq '' or $have_keep = 1; 
	$parm{'drop'} eq '' or $have_drop = 1; 
	my %good_flag; 
	if ( $have_keep == 0 and $have_drop == 0 ) {
		for (sort {$a<=>$b} keys %flag) {
			$good_flag{$_} = 1; 
		}
		return \%good_flag; 
	}
	FLAG_NUM: 
	for my $flag_num ( sort {$a<=>$b} keys %flag ) {
		my @ta = @{$flag{$flag_num}}; 
		my ($drop_flag, $keep_flag) = (0, 0); 
		if ( $have_drop ) {
			CHK_TTL_DROP: 
			for my $exp1 ( split(/;/, $parm{'drop'}) ) {
				my $should_drop = 1; 
				CHK_EACH_DROP: 
				for my $exp2 ( split(/,/, $exp1) ) {
					$exp2 =~ m!^(\d+)=([01])$! or die "Failed to parse |$exp2|\n"; 
					my ($coln, $colv)  = ($1, $2); 
					$ta[ $coln ] != $colv and do { $should_drop = 0; last CHK_EACH_DROP; }; 
				}
				$should_drop == 1 and do { $drop_flag = 1; last CHK_TTL_DROP; }; 
			}
		}
		if ( $have_keep ) {
			CHK_TTL_KEEP: 
			for my $exp1 ( split(/;/, $parm{'keep'}) ) {
				my $should_keep = 1; 
				CHK_EACH_KEEP: 
				for my $exp2 ( split(/,/, $exp1) ) {
					$exp2 =~ m!^(\d+)=([01])$! or die "Failed to parse |$exp2|\n"; 
					my ($coln, $colv) = ($1, $2); 
					$ta[ $coln ] != $colv and do { $should_keep = 0; last CHK_EACH_KEEP; }; 
				}
				$should_keep == 1 and do { $keep_flag = 1; last CHK_TTL_KEEP; }; 
			}
		}
		$have_drop and $drop_flag == 1 and next FLAG_NUM; 
		$have_keep and $keep_flag == 0 and next FLAG_NUM; 
		$good_flag{ $flag_num } = 1; 
	}
	return \%good_flag; 
}# mk_flag()

# Make FLAG list
sub _sam_flag_href {
	my %flag; 
	for ( 0 .. 4095 ) {
		my $binum = unpack( "B32", pack("N", $_) ); 
		$flag{$_} = [ reverse( split(//, sprintf("%012d", $binum) ) ) ]; 
	}
	return \%flag; 
}# _sam_flag_href ()


=head1 sam_flag_infor( $sam_flag )

Function  : 

Return    : (\@flag_infor)

@ [ 0 - 10 ] = [ 0/1, 'explanation' ]; 

$infor_flag[ 0 ] = [qw/0     0x0001  p/, "the read is paired in sequencing"];

$infor_flag[ 1 ] = [qw/1     0x0002  P/, "the read is mapped in a proper pair"];

$infor_flag[ 2 ] = [qw/2     0x0004  u/, "the query sequence itself is unmapped"];

$infor_flag[ 3 ] = [qw/3     0x0008  U/, "the mate is unmapped"];

$infor_flag[ 4 ] = [qw/4     0x0010  r/, "strand of the query (1 for reverse)"];

$infor_flag[ 5 ] = [qw/5     0x0020  R/, "strand of the mate (1 for reverse)"];

$infor_flag[ 6 ] = [qw/6     0x0040  1/, "the read is the first read in a pair"];

$infor_flag[ 7 ] = [qw/7     0x0080  2/, "the read is the second read in a pair"];

$infor_flag[ 8 ] = [qw/8     0x0100  s/, "the alignment is not primary"];

$infor_flag[ 9 ] = [qw/9     0x0200  f/, "the read fails platform/vendor quality checks"];

$infor_flag[ 10] = [qw/10    0x0400  d/, "the read is either a PCR or an optical duplicate"];

$infor_flag{ 11} = [qw/11    0x0800  d/, "the alignment is supplementary alignment"]; 

=cut
sub sam_flag_infor {
	my $flag_num = shift; 
	$flag_num =~ m/^\d+$/ or &stopErr("[Err] Unknown flag_num [$flag_num] in &sam_flag_infor()\n"); 
	$flag_num <= 4095 or &stopErr("sam_flag_number should not be larger than 4095.\n"); 
	my $binum = unpack("B32", pack("N", $flag_num)); 
	my @flag_arr = reverse( split(//, sprintf("%012d", $binum) ) ); 
	my @flag_back; 
	for (my $i=0; $i<@flag_arr; $i++) {
		$flag_back[$i] = [ $flag_arr[$i], $infor_flag{$i}[3] ]; 
	}
	return ( \@flag_back ); 
}

=head1 cigar_str2array ( $CigarString_inSam )

Return      : (\@cigar_num_array)
 @cigar_num_array = ( [cigar_len, cigar_tag], [cigar_len, cigar_tag], ... )

=cut
sub cigar_str2array {
	my $str = shift; 
	my @back; 
	$str eq '*' and do { $str = ''; push(@back, [-1, 'N']);  }; 
	while ($str =~ s!^(\d+)([MIDNSHP=X])!!) {
		push(@back, [$1, $2]); 
	}
	$str eq '' or &stopErr("[Err] Left unknown Cigar string [$str]\n"); 
	return \@back; 
}# sub cigar_str2array () 

=head1 cigar_array2str ( \@cigar_num_array )

Return      :  ($CigarString_inSam)

=cut
sub cigar_array2str {
	my $back = ''; 
	defined $_[0] and @{$_[0]} > 0 and $back = join('', map { "$_->[0]$_->[1]" } @{$_[0]}); 
	return $back; 
}# cigar_array2str ()


=head1 cigar_array2len ( \@cigar_array )

Return      : ( $rd_len, $ref_span_len )

=cut 
sub cigar_array2len {
	my ($rd_len, $ref_len) = (0, 0); 
	for ( @{$_[0]} ) {
		$_->[1] =~ m/^[MIS=X]$/i and $rd_len += $_->[0]; 
		$_->[1] =~ m/^[MDNP=X]$/i and $ref_len += $_->[0]; 
	}
	return ($rd_len, $ref_len); 
}# sub cigar_array2len() 


=head1 trim_pos ( posi_trimmed_before, end_posi_for_trim, add_len_for_trim )

Return      : ( posi_trimmed, left_len_for_trim ) 

=cut
sub trim_pos ($$$) {
	# ( posi_trimmed_before, end_posi_for_trim, add_len_for_trim )
	#   0                    1                  2
	if      ( $_[1]-$_[0] < 1  ) {
		# 1 to end_posi have already been trimmed, so there is no trimming. 
		return($_[0], $_[2]); 
	} elsif ( $_[1]-$_[0] >= $_[2] ) {
		# All of the add_len_for_trim should be trimmed. 
		return($_[0]+$_[2], 0); 
	} elsif ( $_[1]-$_[0] <  $_[2] ) {
		# The heading ($_[1]-$_[0]) positions should be trimmed, and left $_[2]-($_[1]-$_[0]) positions. 
		return($_[1], $_[2]+$_[0]-$_[1]); 
	} else {
		&stopErr("[Err] Why here [@_]\n"); 
	}
}# sub trim_pos () 

=head1 trim_cigar_arr (\@cigar_array, $trim_rd_len)

@cigar_array will be trimmed. 

Return      : ($trimmed_rd_len, $trimmed_ref_len)

=cut
sub trim_cigar_arr ($$) {
	# ( [ [len_1, tag_1], [len_2, tag_2], ... ] , $rd_len_to_trim )
	my @new_cigar_arr; 
	my $trimmed_rd_len = 0; 
	my $trimmed_ref_len = 0; 
	for (@{$_[0]}) {
		if ( $trimmed_rd_len >= $_[1] ) {
			push(@new_cigar_arr, $_); 
			next; 
		}
		if ( $_->[0] == -1 and $_->[1] eq 'N' ) {
			; 
		} elsif ( $_->[1] =~ m/^[MP=X]$/i ) {
			my ($trimmed_RdP, $rest_RdBp) = &trim_pos( $trimmed_rd_len, $_[1], $_->[0] ); 
			$trimmed_ref_len += ($trimmed_RdP - $trimmed_rd_len); 
			$trimmed_rd_len = $trimmed_RdP; 
			$rest_RdBp > 0 and push(@new_cigar_arr, [ $rest_RdBp, $_->[1] ]); 
		} elsif ( $_->[1] =~ m/^I$/i ) {
			$trimmed_rd_len += $_->[0]; 
		} elsif ( $_->[1] =~ m/^S$/i) {
			my ($trimmed_RdP, $rest_RdBp) = &trim_pos( $trimmed_rd_len, $_[1], $_->[0] );
			$trimmed_rd_len = $trimmed_RdP; 
			$rest_RdBp > 0 and push(@new_cigar_arr, [ $rest_RdBp, $_->[1] ]); 
		} elsif ( $_->[1] =~ m/^D$/i) {
			$trimmed_ref_len += $_->[0]; 
		} else {
			; 
		}
	}
	@{$_[0]} = @new_cigar_arr; 
	return ( $trimmed_rd_len, $trimmed_ref_len ); 
}# sub trim_cigar_arr ()

=head1 trim_cigar_str_bothEnd ( $cigar_string, $length_to_trim_from_end ) 

Return      : ( \@trimmed_cigar_array, $left_rm_RdBpN, $right_rm_RdBpN, $left_rm_RefBpN, $right_rm_RefBpN )

=cut 
sub trim_cigar_str_bothEnd {
	# ( $cigar_string , $length_to_trim_from_end )
	my @cigar_arr = @{ &cigar_str2array( $_[0] ) }; 
	# my ($raw_rd_len, $raw_refSpan_len) = &cigar_array2len(\@cigar_arr); 
	
	# Check from the left end. 
	my ( $left_rm_RdBpN,  $left_rm_RefBpN  ) = &trim_cigar_arr( \@cigar_arr, $_[1] ); 
	# Check from the right end. 
	@cigar_arr = reverse(@cigar_arr); 
	my ( $right_rm_RdBpN, $right_rm_RefBpN ) = &trim_cigar_arr( \@cigar_arr, $_[1] ); 
	@cigar_arr = reverse(@cigar_arr); 

	return(\@cigar_arr, $left_rm_RdBpN, $right_rm_RdBpN, $left_rm_RefBpN, $right_rm_RefBpN); 
}

=head1 parseCigar( $CigarString_inSam )

Return      : (\%type2Number)
  {qw/Dlen Elen Hlen Ilen Mlen Nlen Plen RdLen Slen SpanRefLen Xlen/} => number; 
  Plen : Read 'Padded SAM' sectioin in manual. standing for '*'. (silent deletion from padded reference)
  Elen : '=', sequence match 
  Xlen : 'X', sequence mismatch 
  Nlen : 'N', this means an intron. If return '-1', it means this read is not aligned, and its cigar is '*', so its 'SpanRefLen' will be '-1' too. 
  Mlen : 'M', alignment match, can be match or mismatch. 
  Ilen : 'I', insertion to the reference. 
  Dlen : 'D', deletion from the reference. 
  Slen : 'S', soft clipping (clipped sequences present in SEQ)
  Hlen : 'H', hard clipping (clipped sequences NOT present in SEQ)
  SpanRefLen : Mlen + Dlen + Nlen + Plen + Elen + Xlen. '-1' for non-aligned read. 
  RdLen      : Mlen + Ilen + Slen + Elen + Xlen
  MatchRdLen : Mlen + Ilen + Elen + Xlen
  MatchRefLen: Mlen + Dlen + Plen + Elen + Xlen # Do not include intron region!!! 

=cut
sub parseCigar {
	my $cigarString = shift;
	my %back;
	for my $tag (qw/M I D N S H P E X/) {
		$back{"${tag}len"} = 0;
	}
	for my $tag (qw/SpanRefLen RdLen/) {
		$back{$tag} = 0;
	}
	if ( $cigarString eq '*' ) {
		$cigarString = ''; 
		$back{'Nlen'} = -1; 
	}
	while ($cigarString =~ s/^(\d+)([MIDNSHP=X])//) {
		# CIGAR String \*|([0-9]+[MIDNSHPX=])+
		my ($len, $tag) = ($1, $2);
		if ( $tag eq 'M' ) {
			$back{'Mlen'} += $len;
		} elsif ( $tag eq 'I' ) {
			$back{'Ilen'} += $len;
		} elsif ( $tag eq 'D' ) {
			$back{'Dlen'} += $len;
		} elsif ( $tag eq 'N' ) {
			$back{'Nlen'} += $len;
		} elsif ( $tag eq 'S' ) {
			$back{'Slen'} += $len;
		} elsif ( $tag eq 'H' ) {
			$back{'Hlen'} += $len;
		} elsif ( $tag eq 'P' ) {
			$back{'Plen'} += $len;
			# Read 'Padded SAM' section in sam format manual. It stands for '*' ;
			&tsmsg("[Wrn] I don't really understand P in cigar string.\n");
		} elsif ( $tag eq '=' ) {
			$back{'Elen'} += $len;
		} elsif ( $tag eq 'X' ) {
			$back{'Xlen'} += $len;
		} else {
			&stopErr("[Err] Unknown Cigar string $len$tag\n");
		}
	}
	$cigarString eq '' or &stopErr("[Err] Left unknown Cigar string $cigarString\n");
	$back{'SpanRefLen'} = $back{'Mlen'} + $back{'Dlen'} + $back{'Nlen'} + $back{'Plen'} + $back{'Elen'} + $back{'Xlen'};
	$back{'RdLen'} = $back{'Mlen'} + $back{'Ilen'} + $back{'Slen'} + $back{'Elen'} + $back{'Xlen'};
	$back{'MatchRdLen'} = $back{'Mlen'} + $back{'Ilen'} + $back{'Elen'} + $back{'Xlen'}; 
	$back{'MatchRefLen'} = $back{'Mlen'} + $back{'Dlen'} + $back{'Plen'} + $back{'Elen'} + $back{'Xlen'}; 
	return \%back;

}# sub parseCigar() 

=head1 sam_line2hash(\@sam_line_array, [@required_infor], $is_ignore_TAGs)

Required    : \@sam_line_array 

Function    : Retrieve information from sam line as many as possible. 
              [@required_infor] check all of 'could be' in Return. 
              read_str mate_str : could be +/-/u , 'u' means not aligned. 
              ins_s ins_e ins_len : will be computed together, but returned separately. These are same to bwa definition. 
               For +/- pairs : ins_s=rd_+_5'1st_mappingBp , ins_e=rd_-_5'1st_mappingBp , ins_len=ins_e-ins_s+1; 
               For -/+ pairs : ins_s=rd_-_5'1st_mappingBp , ins_e=rd_+_5'1st_mappingBp , ins_len=ins_e-ins_s-1; 
               For +/+ pairs : ins_s=rd_+Curr_5'1st_mappingBp , ins_e=rd_+Mate_5'1st_mappingBp , ins_len=ins_e-ins_s; 
               For -/- pairs : ins_s=rd_-Curr_3'1st_mappingBp , ins_e=rd_-Mate_3'1st_mappingBp , ins_len=ins_e-ins_s; 
               If any end is 'u', all ins_values are 'u' ; 
              
              Set $is_ignore_TAGs as 1 to ignore TAG informations from column_12 for speed. But This requires @required_infor as empty. 


Return      : (\%infor_hash)
 In this %infor_hash, there are keys as 
  qw/qname flag rname pos mapq cigar rnext pnext temlen seq qual/ ; 
  could be [qw/read_len nm_ratio read_str mate_str ins_s ins_e ins_len cigar_href/] ; 
  could be [qw/XA_aref XA_minNM is_uniqBest NM XT XA/] 
    In default : 
    cigar_href = &parseCigar( cigarString ) 
    read_len   = cigar_href->{'RdLen'}
    NM         = 0 ; 
    XT         = 'U' ; 
    XA         = '' ; 
    nm_ratio   = NM / read_len
    read_str|mate_str = +|-|u
    ins_s|ins_e|ins_len = 'u' 
    XA_aref    = [] ; Format [ [$chr_1, $pos_1, $cigar_1, $nm_1], ... ] 
    XA_minNM   = 999999 ; 
    is_uniqBest = 1|0 ; Similar to not_uniqBest(), rules: Not contain 'XT:A:R' + XA:*NM > NM 
    hDiff_Pair = 1|0 ; Set as 1 if both ends mapped, and one is forward and the other reverse . 
=cut

sub sam_line2hash {
	my ($line_aref, $require_aref, $only_basic) = @_; 
	$require_aref //= []; 
	$only_basic //= 0; 
	my %back; 
	@back{qw/qname flag rname pos mapq cigar rnext pnext temlen seq qual/} = @{$line_aref}[0 .. 10]; 
	if (@$line_aref > 11 and !($only_basic and @$require_aref == 0)) {
		for my $tb (@{$line_aref}[11 .. $#$line_aref]) {
			$tb =~ m/^([^:\s]+):([^:\s]+):(.*)$/ or do { &tsmsg("[Wrn] Skip bad tag:value pair [$tb]\n"); }; 
			my ($tag, $type, $value) = ($1, $2, $3); 
			defined $back{$tag} and do { &tsmsg("[Wrn] Skip repeated TAG [$tag]\n"); }; 
			$back{$tag} = $value; 
		}
	}
	my %inner = %back; 
	for my $want (@$require_aref) {
		&_calc_sam_required( \%inner, $want ); 
		$back{$want} = $inner{$want}; 
	}
	return (\%back); 
}# sam_line2hash()

=head1 sam_hash_addKey(\%out_of_sam_line2hash, [@required_infor])

Function   : Add keys in @required_infor to %out_of_sam_line2hash ; 
             %out_of_sam_line2hash is output of sam_line2hash(); 
             And this sub-routine should follow sam_line2hash(); 
=cut
sub sam_hash_addKey {
	my ($line_href, $require_aref) = @_; 
	defined $require_aref or return; 
	my %inner = %$line_href; 
	for my $want (@$require_aref) {
		&_calc_sam_required( \%inner, $want ); 
		$line_href->{$want} = $inner{$want}; 
	}
	return ; 
}#sam_hash_addKey()


# _calc_sam_required( \%basic_hash_from_sam_line2hash, $key_to_set )
# This will change %basic_hash_from_sam_line2hash , and add key($key_to_set) to %basic_hash_from_sam_line2hash; 
# Accept key_need : 
#   cigar_href read_len nm_ratio read_str mate_str ins_s ins_e ins_len 
#   XA_aref XA_minNM is_uniqBest NM XT XA 
# Default : 
#   NM //= 0 ; XT //= 'U' ; XA //= ''; 
#   XA_aref //= [] ; Format [ [$chr_1, $pos_1, $cigar_1, $nm_1], ... ]
#   XA_minNM //= 999999 ; 
#   read_str/mate_str = +|-|u
#   ins_s/ins_e/ins_len = \d+|u
#   is_uniqBest = 1|0
#   hDiff_Pair = 1|0
sub _calc_sam_required {
	my ($inner_href, $key_need) = @_; 
	defined $inner_href->{$key_need} and return; 
	if ( $key_need eq 'cigar_href' ) {
		$inner_href->{$key_need} = &parseCigar( $inner_href->{'cigar'} ); 
	} elsif ( $key_need eq 'read_len' ) {
		&_calc_sam_required($inner_href, 'cigar_href'); 
		$inner_href->{$key_need} = $inner_href->{'cigar_href'}{'RdLen'}; 
	} elsif ( $key_need eq 'nm_ratio' ) {
		&_calc_sam_required($inner_href, 'read_len'); 
		&_calc_sam_required($inner_href, 'NM'); 
		$inner_href->{$key_need} = $inner_href->{'NM'} / $inner_href->{'read_len'}; 
	} elsif ( $key_need eq 'read_str' ) {
		&_chk_sam_flag( 'base_2', $inner_href->{'flag'} ) and do { $inner_href->{'read_str'} = 'u'; return; }; 
		$inner_href->{$key_need} = ( &_chk_sam_flag( 'base_4', $inner_href->{'flag'} ) ) ? '-' : '+' ; 
	} elsif ( $key_need eq 'mate_str' ) {
		&_chk_sam_flag( 'base_3', $inner_href->{'flag'} ) and do { $inner_href->{'mate_str'} = 'u'; return; }; 
		$inner_href->{$key_need} =  ( &_chk_sam_flag( 'base_5', $inner_href->{'flag'} ) ) ? '-' : '+' ; 
	} elsif ( $key_need eq 'ins_s' or $key_need eq 'ins_e' or $key_need eq 'ins_len' ) {
		( &_chk_sam_flag( 'pair_aligned', $inner_href->{'flag'} ) and $inner_href->{'rnext'} eq '=' and $inner_href->{'pnext'} ne '*' ) or do { $inner_href->{'ins_s'}='u'; $inner_href->{'ins_e'}='u'; $inner_href->{'ins_len'}='u'; return; }; 
		&_calc_sam_required($inner_href, 'read_str'); 
		&_calc_sam_required($inner_href, 'mate_str'); 
		if ( $inner_href->{'read_str'} eq '+' ) {
			if ( $inner_href->{'mate_str'} eq '-' ) {
				$inner_href->{'ins_s'} = $inner_href->{'pos'}; 
				$inner_href->{'ins_e'} = $inner_href->{'pos'}+$inner_href->{'temlen'}-1; 
				$inner_href->{'ins_len'} = $inner_href->{'temlen'}; 
				return; 
			} else {
				$inner_href->{'ins_s'} = $inner_href->{'pos'}; 
				$inner_href->{'ins_e'} = $inner_href->{'pnext'}; 
				$inner_href->{'ins_len'} = $inner_href->{'ins_e'}-$inner_href->{'ins_s'}; 
				$inner_href->{'ins_len'} == $inner_href->{'temlen'} or &tsmsg("[Wrn] INS_LEN diff for [$inner_href->{'qname'}] in [$inner_href->{'rname'}]\n"); 
			}
		} else {
			if ( $inner_href->{'mate_str'} eq '+' ) {
				$inner_href->{'ins_e'} = $inner_href->{'pnext'}; 
				$inner_href->{'ins_s'} = $inner_href->{'ins_e'} - $inner_href->{'temlen'} + 1; 
				$inner_href->{'ins_len'} = $inner_href->{'temlen'}; 
			} else {
				$inner_href->{'ins_s'} = $inner_href->{'pos'}; 
				$inner_href->{'ins_e'} = $inner_href->{'pnext'}; 
				$inner_href->{'ins_len'} = $inner_href->{'ins_e'}-$inner_href->{'ins_s'}; 
				$inner_href->{'ins_len'} == $inner_href->{'temlen'} or &tsmsg("[Wrn] INS_LEN diff for [$inner_href->{'qname'}] in [$inner_href->{'rname'}]\n"); 
			}
		}
	} elsif ( $key_need eq 'XA_aref' ) {
		my @t_xa_aref; 
		&_calc_sam_required($inner_href, 'XA'); 
		for my $ts0 (split(/;/, $inner_href->{'XA'})) {
			$ts0 =~ m/^\s*$/ and next; 
			$ts0 =~ m/^(\S+),([+-]?\d+),([\d\w]+),(\d+)$/ or &stopErr("[Err] Failed for XA:Z : [$ts0]\n"); 
			my ($chr, $pos, $cigar, $nm_1) = ($1, $2, $3, $4); 
			push(@t_xa_aref, [$chr, $pos, $cigar, $nm_1]); 
		}
		$inner_href->{$key_need} = \@t_xa_aref; 
	} elsif ( $key_need eq 'XA_minNM' ) {
		&_calc_sam_required($inner_href, 'XA_aref'); 
		my $xa = 999999; 
		for my $ar (@{$inner_href->{'XA_aref'}}) {
			$ar->[3] < $xa and $xa = $ar->[3]; 
		}
		$inner_href->{$key_need} = $xa; 
	} elsif ( $key_need eq 'is_uniqBest' ) {
		&_calc_sam_required($inner_href, 'XT'); 
		$inner_href->{'XT'} eq 'R' and do { $inner_href->{$key_need} = 0; return; }; 
		&_calc_sam_required($inner_href, 'NM'); 
		&_calc_sam_required($inner_href, 'XA_minNM'); 
		$inner_href->{'NM'} >= $inner_href->{'XA_minNM'} and do { $inner_href->{$key_need} = 0; return; }; 
		$inner_href->{$key_need} = 1; 
	} elsif ( $key_need eq 'hDiff_Pair' ) {
		$inner_href->{$key_need} = ( &_chk_sam_flag( 'hDiff_Pair', $inner_href->{'flag'} ) ) ? 1 : 0 ; 
	} elsif ( $key_need eq 'NM' ) {
		$inner_href->{$key_need} = 0; 
	} elsif ( $key_need eq 'XT' ) {
		$inner_href->{$key_need} = 'U'; 
	} elsif ( $key_need eq 'XA' ) {
		$inner_href->{$key_need} = ''; 
	} else {
		&stopErr("[Err] Why here? key_need=[$key_need]\n"); 
	} 
	return; 
} # _calc_sam_required()

=head1 print_sam_lines( \@rdID_wiLine, \%rdID_toDrop, $outfile_handle ) 

Function    : Print out sam_lines in @rdID_wiLine after removing by %rdID_toDrop with file handle $outfile_handle ; 

Input       : 
  @rdID_wiLine : ( [$rdID, $line_without_return], [$rdID, $line_without_return], ... )
  %rdID_toDrop : ( $rdID => 1 , $rdID => 1 ) ; Any rdID defined in this hash won't be output. 
  $outfile_handle : Default \*STDOUT 

Return      : $out_line_number 

=cut
sub print_sam_lines {
	my ($line_aref, $toRM_href, $out_fh) = @_; 
	$toRM_href //= {}; 
	$out_fh //= \*STDOUT; 
	my $line_num = 0; 
	for my $lar (@$line_aref) {
		defined $toRM_href->{$lar->[0]} and next; 
		print {$out_fh} $lar->[1] . "\n"; 
		$line_num ++; 
	}
	return $line_num; 
}#print_sam_lines() 


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
	unless ( ref($seqA) eq 'SCALAR' or ref($seqA) eq '' ) {
		ref($seqA) eq 'SeqAlnSunhh' or &stopErr("[Err] olap_e2e_A2B 1st-input should be a scalar string standing for sequence A.\n", "Now ref_seqA is [", ref($seqA), "]\n", "seqA=$seqA\n"); 
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

