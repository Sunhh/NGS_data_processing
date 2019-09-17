#!/usr/bin/perl
# Honghe Sun : sunhonghe@nercv.org
# Updated    : 20190916
use strict; 
use warnings; 
use LogInforSunhh; 
use fastaSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"task:s",        # "pcr" supported; 
	"in_loc:s",     # location table : locID \\t seqID \\t taget_start \\t target_end \\t temp_start \\t temp_end 
	                 #   In default, [temp_start, temp_end] = [target_start - 1000, target_end + 1000]
	"in_seq:s",      # fasta file with seqID; 
	"in_p3conf:s",   # optional; directly used in primer3 input file, with necessary edit by sconf_* parameters. 
	"sconf_primer_size:s",   # opt:min:max; 
	"sconf_product_size:s@", # length ranges 
	"sconf_return_num:i",      # number of pairs returned; 
	"exe_primer3:s",         # '/home/Sunhh/bin/primer3_core'
	"default_version:i",     # 2 
	"sconf_th_path:s",             # /Data/Sunhh/src/genetics/primer3/primer3-2.4.0/src/primer3_config/ 
	"help!", 
); 

my $h_txt = <<HH; 
######################################################################
# perl $0 -in_seq in.fas  -in_loc tgt.loc  > out
#
# -in_seq               fasta file. 
# -in_loc               locID \\t seqID \\t taget_start \\t target_end \\t temp_start \\t temp_end
# -in_p3conf            filename. similar to primer3 input file. 
#
# -sconf_primer_size    [opt:min:max]
# -sconf_product_size   [num-num] could be multiple times. 
# -sconf_return_num     [num] for PRIMER_NUM_RETURN
# -in_p3conf            filename. Additional parameters for primer3; 
#
# -exe_primer3          path to primer3_core
# -default_version      Default 2 for the latest version; 1 for old; 
# -sconf_th_path        path to thermodynamic parameters; 
# 
# -help 
######################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($h_txt); 
defined $opts{'in_seq'} or &LogInforSunhh::usage($h_txt); 
defined $opts{'in_loc'} or &LogInforSunhh::usage($h_txt); 

my $fs_obj = fastaSunhh->new(); 

my %gg; 
### Set up %gg from %opts; 
for my $k1 (sort keys %opts) {
	$gg{$k1} = $opts{$k1}; 
}
### Set up %gg from in_p3conf file.
if (defined $gg{'in_p3conf'}) {
	my $ifh=&openFH($gg{'in_p3conf'}, '<'); 
	### Multiple tags not allowed! 
	while (<$ifh>) {
		chomp; 
		m!^\s*=\s*$! and next; 
		m!^\s*$! and next; 
		m!^(\S+)=(.+)$! or &stopErr("[Err] bad line: $_\n"); 
		my ($tag, $key) = ($1, $2); 
		$gg{'inconf'}{$tag} = $key; 
	}
	close($ifh); 
}
### Update 'inconf' TAGs assigned by %opts sconf; 
if (defined $gg{'sconf_primer_size'}) {
	$gg{'sconf_primer_size'} =~ m!^([\d\.]+):([\d\.]+):([\d\.]+)$! or &stopErr("[Err] Bad value [$gg{'sconf_primer_size'}] for sconf_primer_size\n "); 
	@{$gg{'inconf'}}{qw/PRIMER_OPT_SIZE PRIMER_MIN_SIZE PRIMER_MAX_SIZE/} = ($1, $2, $3); 
}
if (defined $gg{'sconf_product_size'}) {
	my @a1; 
	for my $v1 (@{$gg{'sconf_product_size'}}) {
		$v1 =~ m!^(\d+)\-(\d+)$! or &stopErr("[Err] Bad value [$v1] for sconf_product_size\n"); 
		push(@a1, $v1); 
	}
	# $gg{'sconf_product_size'} = join(" ", @a1); 
	$gg{'inconf'}{'PRIMER_PRODUCT_SIZE_RANGE'} = join(" ", @a1); 
}
if (defined $gg{'sconf_return_num'}) {
	$gg{'inconf'}{'PRIMER_NUM_RETURN'} = $gg{'sconf_return_num'}; 
}
if (defined $gg{'sconf_th_path'}) {
	$gg{'inconf'}{'PRIMER_THERMODYNAMIC_PARAMETERS_PATH'} = $gg{'sconf_th_path'}; 
}
### Load fasta templates and target loci; 
$gg{'dbFa'} = $fs_obj->save_seq_to_hash( 'faFile' => $gg{'in_seq'} ); 
for my $tk (keys %{$gg{'dbFa'}}) { $gg{'dbFa'}{$tk}{'seq'} =~ s!\s+!!g; $gg{'dbFa'}{$tk}{'len'} = length($gg{'dbFa'}{$tk}{'seq'}); }
for my $tl (&fileSunhh::load_tabFile($gg{'in_loc'})) {
	$tl->[0] eq 'locID' and next; 
	(defined $tl->[4] and $tl->[4] !~ m!^(\s*|NA)$!i) or $tl->[4] = $tl->[2] - 1000; # temp_start
	(defined $tl->[5] and $tl->[5] !~ m!^(\s*|NA)$!i) or $tl->[5] = $tl->[3] + 1000; # temp_end
	# Check if range exceeds the length of template sequence. 
	$tl->[4] < 1 and $tl->[4] = 1; 
	$tl->[5] > $gg{'dbFa'}{$tl->[1]}{'len'} and $tl->[5] = $gg{'dbFa'}{$tl->[1]}{'len'}; 
	$tl->[4] <= $tl->[2] or &stopErr("[Err] Bad temp_start [$tl->[4]] in line: @$tl\n"); 
	$tl->[5] >= $tl->[3] or &stopErr("[Err] Bad temp_end [$tl->[5]] in line: @$tl\n"); 
	$tl->[5] > $tl->[4] or &stopErr("[Err] temp_start [$tl->[4]] should be smaller than temp_end [$tl->[5]] in line: @$tl\n"); 
	my $tseq = substr($gg{'dbFa'}{$tl->[1]}{'seq'}, $tl->[4]-1, $tl->[5]-$tl->[4]+1); 
	my $tgt_s = $tl->[2]-$tl->[4]+1; 
	my $tgt_e = $tl->[3]-$tl->[4]+1; 
	my $tgt_l = $tgt_e - $tgt_s + 1; 
	push(@{$gg{'tgt'}{'loc'}}, [@{$tl}[0..5], $tseq, $tgt_s, $tgt_e, $tgt_l]); # locID  seqID  taget_start  target_end  temp_start  temp_end  temp_seq tgtS tgtE tgtLen
}

### Complete required %gg; 
$gg{'exe_primer3'} //= '/Data/Sunhh/src/genetics/primer3/primer3-2.4.0/src/primer3_core'; 
$gg{'default_version'} //= 2; 
$gg{'inconf'}{'PRIMER_THERMODYNAMIC_PARAMETERS_PATH'} //= '/Data/Sunhh/src/genetics/primer3/primer3-2.4.0/src/primer3_config/'; 
# $gg{'inconf'}{'PRIMER_FIRST_BASE_INDEX'} = 1; 

### Run primer3 one by one and output results; 
print STDOUT join("\t", qw/locID markerID l_seq r_seq l_tm r_tm l_penalty r_penalty l_5pPos l_3pPos r_3pPos r_5pPos l_gc r_gc product_size product_seq/)."\n"; 
my $wdir = &fileSunhh::new_tmp_dir('create' => 1); 
for my $t1 (@{$gg{'tgt'}{'loc'}}) {
	my $ofh1 = &openFH("$wdir/input", '>'); 
	print {$ofh1} "SEQUENCE_ID=$t1->[0]\n"; 
	print {$ofh1} "SEQUENCE_TEMPLATE=$t1->[6]\n"; 
	print {$ofh1} "PRIMER_FIRST_BASE_INDEX=1\n"; 
	print {$ofh1} "SEQUENCE_TARGET=$t1->[7],$t1->[9]\n"; 
	for my $t2 (sort keys %{$gg{'inconf'}}) {
		$t2 eq 'SEQUENCE_ID' and next; 
		$t2 eq 'SEQUENCE_TEMPLATE' and next; 
		$t2 eq 'SEQUENCE_TARGET' and next; 
		$t2 eq 'PRIMER_FIRST_BASE_INDEX' and next; 
		print {$ofh1} "$t2=$gg{'inconf'}{$t2}\n"; 
	}
	print {$ofh1} "=\n"; 
	close($ofh1); 
	&run_cmd("$gg{'exe_primer3'} --default_version=$gg{'default_version'} --strict_tags --output=$wdir/p3_out --error=$wdir/p3_err $wdir/input"); 
	# Load output; 
	my $tfh1 = &openFH("$wdir/p3_out", '<'); 
	my @p3out = ({}); 
	while (<$tfh1>) {
		chomp; 
		m!^\=$! and do { push(@p3out, {}); next; }; 
		m!^(\S+?)=(.+)$! or &stopErr("[Err] in file $wdir/p3_out bad line: $_\n"); 
		$p3out[-1]{$1} = $2; 
	}
	close($tfh1); 
	# Output output table with in_loc information; 
	for my $t2 (@p3out) {
		scalar(keys %$t2) == 0 and next; 
		RET_PRIMER: 
		for (my $pv=0; 1; $pv++) {
			my (%l_p, %r_p, %i_p); # left, right, internal primer/oligo; 
			# For paired primers only; 
			my %t3_seq = %{ &fill_primer_v($pv, "SEQUENCE", $t2) }; 
			defined $t3_seq{'LEFT'}  or last RET_PRIMER; 
			defined $t3_seq{'RIGHT'} or last RET_PRIMER; 
			my %t3_tm   = %{ &fill_primer_v($pv, "TM", $t2) }; 
			my %t3_loc  = %{ &fill_primer_v($pv, "", $t2) }; 
			my %t3_pen  = %{ &fill_primer_v($pv, "PENALTY", $t2) }; 
			my %t3_gc   = %{ &fill_primer_v($pv, "GC_PERCENT", $t2) }; 
			my %t3_sany = %{ &fill_primer_v($pv, "SELF_ANY", $t2) }; 
			my %t3_send = %{ &fill_primer_v($pv, "SELF_END", $t2) }; 
			my %t3_estab = %{ &fill_primer_v($pv, "END_STABILITY", $t2) }; 
			my %t3_cany  = %{ &fill_primer_v($pv, "COMPL_ANY", $t2) }; # pair
			my %t3_cend  = %{ &fill_primer_v($pv, "COMPL_END", $t2) }; # pair 
			my %t3_psize = %{ &fill_primer_v($pv, "PRODUCT_SIZE", $t2) }; # pair 

			my ($l_5p, $l_3p, $r_5p, $r_3p); 
			$t3_loc{'LEFT'} =~ m!^(\d+),(\d+)$! or &stopErr("[Err] bad LEFT loc $t3_loc{'LEFT'}\n"); 
			my @lsl = ($1, $2); 
			$l_5p = $lsl[0] + $t1->[4] - 1; 
			$l_3p = $l_5p + $lsl[1] - 1; 
			$t3_loc{'RIGHT'} =~ m!^(\d+),(\d+)$! or &stopErr("[Err] bad RIGHT loc $t3_loc{'RIGHT'}\n"); 
			my @rsl = ($1, $2); 
			$r_5p = $rsl[0] + $t1->[4] - 1; 
			$r_3p = $r_5p - $rsl[1] + 1; 
			my $product_seq = substr($t2->{'SEQUENCE_TEMPLATE'}, $lsl[0]-1, $rsl[0]-$lsl[0]+1); 

			print STDOUT join("\t", 
				$t1->[0],         # locID
				"$t1->[0]_${pv}", # locID_primerNum
				$t3_seq{'LEFT'},  # left seq
				$t3_seq{'RIGHT'}, # right seq
				$t3_tm{'LEFT'},   # left TM
				$t3_tm{'RIGHT'},  # right TM
				$t3_pen{'LEFT'},  # left penalty
				$t3_pen{'RIGHT'}, # right penalty
				$l_5p, # left 5' in in_fa
				$l_3p, # left 3' in in_fa
				$r_3p, # right 3' in in_fa
				$r_5p, # right 5' in in_fa
				$t3_gc{'LEFT'}, # left GC
				$t3_gc{'RIGHT'}, # right GC
				$t3_psize{'PAIR'}, # product size
				$product_seq, # product sequence
			)."\n"; 
		}
	}
}

&fileSunhh::_rmtree("$wdir"); 


sub fill_primer_v{
	my ($pv, $suff, $th) = @_; 
	my %back; 
	for my $c (qw/LEFT RIGHT INTERNAL PAIR/) {
		my $t = "PRIMER_${c}_${pv}_$suff"; 
		if ($suff eq "") {
			$t = "PRIMER_${c}_${pv}"; 
		}
		defined $th->{$t} or next; 
		defined $back{$c} and &stopErr("[Err] Repeat for back{$c}\n"); 
		$back{$c} = $th->{$t}; 
	}
	return(\%back); 
}# fill_primer_v 

sub run_cmd {
	&exeCmd_1cmd($_[0]) and &stopErr("[Err] Failed at CMD: $_[0]\n"); 
}

