#!/usr/bin/perl
# 2018-10-25 : Retrieve read alignments according to src2tgt_rdList from _comb.bam; 
#   Input   : wrk_dir/pref_comb.bam
#   Input   : wrk_dir/pref.src2tgt_rdList
#   Produce : wrk_dir/pref.src2tgt_rd.bam
use strict; 
use warnings; 
use LogInforSunhh; 
use SeqAlnSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"pref:s",    # Project prefix; 
		"inBam:s", # Input bam file, which comes from step1.pl 
		"inRdList:s", 
		"outBam:s", 
	"wrk_dir:s",        # './', Assign a directory to store all resulting bam files, instead of in current folder. 
	
	"exe_samtools:s", 
); 

my %flag_UN = %{ &SeqAlnSunhh::mk_flag( 'keep' => '2=1' ) }; 
my %gg; 
&setGlob(); 
&applyOpt(); 
&get_bam(); 
sub get_bam {
	-e $gg{'inRdList'} or &stopErr("[Err] src2tgt_rdList -inRdList [$gg{'inRdList'}] not exists.\n"); 
	-e $gg{'inBam'} or &stopErr("[Err] Failed to find input _comb.bam file [$gg{'inBam'}]\n"); 
	# my $outBam = "$gg{'wrk_dir'}/$gg{'pref'}.src2tgt_rd.bam"; 

	my $ifh = &openFH($gg{'inRdList'}, '<'); 
	my %maxN; 
	# readID src_misN tgt_misN
	while (<$ifh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		$ta[1] //= -1; 
		$ta[1] eq 'src_misN' and next; 
		$maxN{$ta[0]} = $ta[1]; 
	}
	close($ifh); 

	my %log_cnt; 
	%log_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>5e6 ); 
	open F,'-|', "$gg{'exe_samtools'} view -h $gg{'inBam'} " or &stopErr("[Err] Failed to read bam file [$gg{'inBam'}]\n"); 
	open O,'|-', "$gg{'exe_samtools'} view -o $gg{'outBam'} - " or &stopErr("[Err] Failed to write bam file [$gg{'outBam'}]\n"); 
	while (<F>) {
		&fileSunhh::log_section($., \%log_cnt) and &tsmsg("[Msg]   Reading $. line from file $gg{'inBam'}\n"); 
		chomp; 
		m!^\@! and do { print O "$_\n"; next; }; 
		my @ta = split(/\t/, $_); 
		defined $maxN{$ta[0]} or next; 
		if ( $maxN{$ta[0]} < 0 ) {
			print O "$_\n"; 
			next; 
		}
		defined $flag_UN{$ta[1]} and next; 
		$ta[5] eq '*' and next; 
		my $statH = &stat_aln(\@ta); 
		$statH->{'whole_mismat'} <= $maxN{$ta[0]} or next; 
		print O "$_\n"; 
	}
	close O; 
	close F; 

	return; 
}# get_bam() 



sub setGlob {
	$gg{'pref'}         = 'out'; 

	$gg{'exe_samtools'} = 'samtools'; 
	$gg{'wrk_dir'}          = './'; 
	$gg{'inRdList'}         = "$gg{'wrk_dir'}/$gg{'pref'}.src2tgt_rdList"; 
	$gg{'inBam'}            = "$gg{'wrk_dir'}/$gg{'pref'}_comb.bam"; 
	$gg{'outBam'}           = "$gg{'wrk_dir'}/$gg{'pref'}.src2tgt_rd.bam"; 

$gg{'help_txt'} = <<"HH"; 
################################################################################
# perl $0   -pref $gg{'pref'} 
#
#   -pref       [$gg{'pref'}] Output prefix 
#     -inBam    [$gg{'inBam'}] Not necessary. This input bam file comes from step1.pl 
#     -inRdList [$gg{'inRdList'}] Not necessary. 
#   -wrk_dir    [$gg{'wrk_dir'}] Directory storing result files. 
#   
#   -outBam      [$gg{'outBam'}]
#
#
#
#   -exe_samtools      [$gg{'exe_samtools'}] 
#
#   Output :  $gg{'outBam'}
################################################################################
HH

	return; 
}# setGlob() 

sub applyOpt {
	$opts{'help'} and &LogInforSunhh::usage($gg{'help_txt'}); 

	# For outer tools
	for my $k (qw/exe_samtools/) {
		defined $opts{$k} and $gg{$k} = $opts{$k}; 
	}

	for my $k (qw/pref wrk_dir inRdList outBam/) {
		defined $opts{$k} and $gg{$k} = $opts{$k}; 
	}
	
	$gg{'inBam'}            = "$gg{'wrk_dir'}/$gg{'pref'}_comb.bam"; 
	defined $opts{'inBam'} and $gg{'inBam'} = $opts{'inBam'}; 
	$gg{'inRdList'}         = "$gg{'wrk_dir'}/$gg{'pref'}.src2tgt_rdList"; 
	defined $opts{'inRdList'} and $gg{'inRdList'} = $opts{'inRdList'}; 
	$gg{'outBam'}           = "$gg{'wrk_dir'}/$gg{'pref'}.src2tgt_rd.bam"; 
	defined $opts{'outBa,'} and $gg{'outBam'} = $opts{'outBam'}; 

	return; 
}# applyOpt() 

sub stat_aln {
	my ($ar) = @_;
	my %backH;
	$backH{'cigarH'} = &SeqAlnSunhh::parseCigar( $ar->[5] );
	$backH{'read_len'} = $backH{'cigarH'}{'RdLen'};
	$backH{'whole_mismat'} = 0;
	for my $tk (qw/Slen Hlen/) {
		defined $backH{'cigarH'}{$tk} and $backH{'whole_mismat'} += $backH{'cigarH'}{$tk};
	}
	for my $tb (@{$ar}[ 11 .. $#$ar ]) {
		$tb =~ m!^NM:i:(\d+)$! or next;
		$backH{'whole_mismat'} += $1;
		last;
	}
	return(\%backH);
}# stat_aln()


