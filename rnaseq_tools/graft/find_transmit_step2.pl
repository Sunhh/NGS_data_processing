#!/usr/bin/perl
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
	"src_fa:s@", # Source genome where transmitted reads come from; 
	"tgt_fa:s@", # Target/Current genome where transmitted reads don't belong to; 
	"wrk_dir:s",        # './', Assign a directory to store all resulting bam files, instead of in current folder. 
	
	"exe_samtools:s", 
); 

my %flag_UN = %{ &SeqAlnSunhh::mk_flag( 'keep' => '2=1' ) }; 
my %gg; 
&setGlob(); 
&applyOpt(); 
&load_allFaID(); # IDs stored in $gg{'refClass'} : {$sourceID} => 'src', {$targetID} => 'tgt'; 
&step2_filter_mismat(); 

sub step2_filter_mismat {
	# %{$gg{'src_fa_rdMisMat'}{$rdID}}; 
	# Produce : 
	#   $gg{'inBam'}.src2tgt_rdList.1 ; 
	#   $gg{'inBam'}.stat    ; 
	
	-e $gg{'inBam'} or &stopErr("[Err] Failed to find input _comb.bam file [$gg{'inBam'}]\n"); 

	my %log_cnt; 
	%log_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>5e6 ); 
	open F,'-|', "$gg{'exe_samtools'} view $gg{'inBam'} " or &stopErr("[Err] Failed to read bam file [$gg{'inBam'}]\n"); 
	my %prevH; 
	while (<F>) {
		&fileSunhh::log_section($., \%log_cnt) and &tsmsg("[Msg]   Reading $. line from file $gg{'inBam'}\n"); 
		chomp; 
		my @ta = split(/\t/, $_); 
		defined $flag_UN{$ta[1]} and next; 
		$ta[5] eq '*' and next; 
		my $statH = &stat_aln(\@ta); 
		my $rdID = "$ta[0]"; 
		defined $gg{'refClass'}{$ta[2]} or do { &tsmsg("[Wrn] Skip unknown refID [$ta[2]]\n"); next; }; 
		my $refClass = $gg{'refClass'}{$ta[2]}; 
		if ( defined $prevH{'rdID'} ) {
			if ( $prevH{'rdID'} eq $rdID ) {
				if ( defined $prevH{'min_mismat'}{$refClass} ) {
					$prevH{'min_mismat'}{$refClass} > $statH->{'whole_mismat'} and $prevH{'min_mismat'}{$refClass} = $statH->{'whole_mismat'}; 
				} else {
					$prevH{'min_mismat'}{$refClass} = $statH->{'whole_mismat'}; 
				}
			} else {
				for my $a0 (qw/tgt src/) {
					$prevH{'min_mismat'}{$a0} //= -1; 
				}
				&process_prevH( \%prevH ); 
				$prevH{'rdID'} = $rdID; 
				$prevH{'min_mismat'}{$refClass} = $statH->{'whole_mismat'}; 
			}
		} else {
			$prevH{'rdID'} = $rdID; 
			$prevH{'min_mismat'}{$refClass} = $statH->{'whole_mismat'}; 
		}
	}
	close F; 
	if ( scalar(keys %prevH) > 0 ) {
		for my $a0 (qw/tgt src/) {
			$prevH{'min_mismat'}{$a0} //= -1; 
		}
		&process_prevH( \%prevH ); 
	}

	# Output potential transmitted reads list. 
	my $ofh_1 = &openFH("$gg{'inBam'}.src2tgt_rdList", '>'); 
	print {$ofh_1} join("\t", qw/readID src_misN tgt_misN/)."\n"; 
	for my $t1 (@{$gg{'list_src2tgt'}}) {
		print {$ofh_1} join("\t", @$t1)."\n"; 
	}
	for my $t0 (qw/mismat_tgtOnly mismat_srcHigh mismat_same mismat_tgtHigh mismat_srcOnly/) {
		$gg{'rdCnt'}{$t0} //= 0; 
	}
	&fileSunhh::write2file("$gg{'inBam'}.stat", "", '>'); 
	&fileSunhh::write2file("$gg{'inBam'}.stat", join("\t", "target_only",   $gg{'rdCnt'}{'mismat_tgtOnly'})."\n", '>>'); 
	&fileSunhh::write2file("$gg{'inBam'}.stat", join("\t", "target_better", $gg{'rdCnt'}{'mismat_srcHigh'})."\n", '>>'); 
	&fileSunhh::write2file("$gg{'inBam'}.stat", join("\t", "same_mismatch", $gg{'rdCnt'}{'mismat_same'})."\n",    '>>'); 
	&fileSunhh::write2file("$gg{'inBam'}.stat", join("\t", "source_better", $gg{'rdCnt'}{'mismat_tgtHigh'})."\n", '>>'); 
	&fileSunhh::write2file("$gg{'inBam'}.stat", join("\t", "source_only",   $gg{'rdCnt'}{'mismat_srcOnly'})."\n", '>>'); 
	
	return; 
}# step2_filter_mismat() 


sub setGlob {
	$gg{'pref'}         = 'out'; 

	$gg{'exe_samtools'} = 'samtools'; 
	$gg{'wrk_dir'}          = './'; 

$gg{'help_txt'} = <<"HH"; 
################################################################################
# perl $0   -pref $gg{'pref'}   -src_fa fa   -tgt_fa fa  
#
#   -pref       [$gg{'pref'}] Output prefix 
#     -inBam    [$gg{'inBam'}] Not necessary. This input bam file comes from step1.pl 
#   -wrk_dir    [$gg{'wrk_dir'}] Directory storing result files. 
#   -src_fa     [faFile] \@ Indexed fasta files for source taxa where transmitted reads come from. 
#   -tgt_fa     [faFile] \@ Indexed fasta files for target taxa where transmitted reads go to. 
#   
#
#
#
#   -exe_samtools      [$gg{'exe_samtools'}] 
################################################################################
HH

	return; 
}# setGlob() 

sub applyOpt {
	$opts{'help'} and &LogInforSunhh::usage($gg{'help_txt'}); 
	for my $k (qw/src_fa tgt_fa/) {
		unless ( defined $opts{$k} ) {
			&tsmsg("[Err] -$k is required.\n\n"); 
			&LogInforSunhh::usage($gg{'help_txt'}); 
		}
	}

	for my $k (qw/pref wrk_dir src_fa tgt_fa/) {
		defined $opts{$k} and $gg{$k} = $opts{$k}; 
	}
	
	$gg{'inBam'}            = "$gg{'wrk_dir'}/$gg{'pref'}_comb.bam"; 
	defined $opts{'inBam'} and $gg{'inBam'} = $opts{'inBam'}; 

	return; 
}# applyOpt() 

sub load_allFaID {
	my %t1; 
	for my $k1 (qw/src tgt/) {
		my $k3 = "${k1}_fa_ID"; 
		for my $fafn (@{$gg{"${k1}_fa"}}) {
			&tsmsg("[Msg] Loading ${k1}_fa file [$fafn]\n"); 
			my $th = &load_faKey( $fafn ); 
			for my $k2 (keys %$th) {
				$t1{$k3}{$k2} = $k1; 
			}
		}
	}

	# Check if there are same ID between source and target. 
	for my $k2 (keys %{$t1{'src_fa_ID'}}) {
		defined $t1{'tgt_fa_ID'}{$k2} and &stopErr("[Err] There is same ID [$k2] between source and target fasta files.\n"); 
	}

	# Store to gg; 
	for my $k1 (qw/src tgt/) {
		my $k3 = "${k1}_fa_ID"; 
		for my $k2 (keys %{$t1{$k3}}) {
			$gg{'refClass'}{$k2} = $k1; 
		}
	}

	return; 
}# load_allFaID() 



sub load_faKey {
	my ($faFn) = @_; 
	my %back; 
	my $fh = &openFH($faFn, '<'); 
	while (<$fh>) {
		m!^\s*\>(\S+)! or next; 
		$back{$1} //= $.; 
	}
	close($fh); 
	return(\%back); 
}# load_faKey

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

sub process_prevH {
	my ($hr) = @_; 
	# There are qw/tgt src rdID/ keys in %$hr; 
	# Save information to 
	#   list_src2tgt => [ rdID, src_mismat, tgt_mismat]
	#   rdCnt=> mismat_srcOnly mismat_tgtOnly mismat_same mismat_srcHigh mismat_tgtHigh => 
	if ($hr->{'min_mismat'}{'tgt'} < 0) {
		# This read may come from source. 
		$gg{'rdCnt'}{'mismat_srcOnly'} ++; 
		push(@{$gg{'list_src2tgt'}}, [ $hr->{'rdID'}, $hr->{'min_mismat'}{'src'}, $hr->{'min_mismat'}{'tgt'} ]); 
	} elsif ($hr->{'min_mismat'}{'src'} < 0) {
		# There is no alignment in source, so I dont want this read. 
		$gg{'rdCnt'}{'mismat_tgtOnly'} ++; 
	} elsif ($hr->{'min_mismat'}{'tgt'} == $hr->{'min_mismat'}{'src'}) {
		$gg{'rdCnt'}{'mismat_same'} ++; 
	} elsif ($hr->{'min_mismat'}{'tgt'} < $hr->{'min_mismat'}{'src'}) {
		$gg{'rdCnt'}{'mismat_srcHigh'} ++; 
	} else {
		# This read may come from source. 
		$gg{'rdCnt'}{'mismat_tgtHigh'} ++; 
		push(@{$gg{'list_src2tgt'}}, [ $hr->{'rdID'}, $hr->{'min_mismat'}{'src'}, $hr->{'min_mismat'}{'tgt'} ]); 
	}
	%{$hr} = (); 
	return; 
}# process_prevH() 


