#!/usr/bin/perl 
# perl /home/Sunhh/tools/github/NGS_data_processing/assemble_tools/format_maf_forMugsy.pl -specs "," S1.maf -out S1.fitMugsy.maf
# perl /home/Sunhh/tools/github/NGS_data_processing/assemble_tools/maf2fasta.pl S1.fitMugsy.maf > S1.fitMugsy.maf.xmfa
# source /data/Sunhh/src/Align/Mugsy/mugsy_x86-64-v1r2.3/mugsyenv.sh ; /data/Sunhh/src/Align/Mugsy/mugsy_x86-64-v1r2.3/mugsyWGA --outfile S1_syn --seq S1_all.fitMugsy.fa --aln S1.fitMugsy.maf.xmfa --distance 2000 --minlength 50 1>stdout.S1_syn 2>stderr.S1_syn
# perl /home/Sunhh/tools/github/NGS_data_processing/assemble_tools/get_paired_maf.pl S1_syn.maf > S1_syn.paired.maf

use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use ReadInAlnSunhh; 
use mathSunhh; 

use Parallel::ForkManager; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"printCmd!",          # Do not execute the final step. 
	"cpuN:i",              # Number of CPUs, default 1 . 
	"seqPerBatch:i",       # Number of sequences in each .maf file, default 200 . Could be slightly larger than this number. 
	"outfile:s",           # Name of out file. Like 'S1_syn' 
	"inmaf:s",             # input file name, like 'S1.fitMugsy.maf'; 
	"infas:s",             # input fasta file, like 'S1_all.fitMugsy.fa'
	
	"pl_maf2fasta:s",      # $HOME/tools/github/NGS_data_processing/assemble_tools/maf2fasta.pl
	"pl_get_paired_maf:s", # $HOME/tools/github/NGS_data_processing/assemble_tools/get_paired_maf.pl
	"pl_format_maf_forMugsy:s", # $HOME/tools/github/NGS_data_processing/assemble_tools/format_maf_forMugsy.pl
	"pl_run_cmd_in_batch:s", # run_cmd_in_batch.pl
	"mugsy_src:s",         # /data/Sunhh/src/Align/Mugsy/mugsy_x86-64-v1r2.3/mugsyenv.sh 
	"mugsy_exe:s",         # /data/Sunhh/src/Align/Mugsy/mugsy_x86-64-v1r2.3/mugsyWGA
	"mugsy_para:s",        # --distance 2000 --minlength 50 
	                       # --outfile S1_syn --seq S1_all.fitMugsy.fa --aln S1.fitMugsy.maf.xmfa
						   # 1>stdout.S1_syn 2>stderr.S1_syn
	"help!", 
); 

sub usage {
	print <<HH; 
################################################################################
# perl $0 -
	"cpuN:i",              # Number of CPUs, default 1 . 
	"seqPerBatch:i",       # Number of sequences in each .maf file, default 1000 . Could be slightly larger than this number. 
	"outfile:s",           # Name of out file. Like 'S1_syn' 
	"inmaf:s",             # input file name, like 'S1.fitMugsy.maf'; 
	
	"pl_maf2fasta:s",      # \$HOME/tools/github/NGS_data_processing/assemble_tools/maf2fasta.pl
	"pl_get_paired_maf:s", # \$HOME/tools/github/NGS_data_processing/assemble_tools/get_paired_maf.pl
	"mugsy_src:s",         # /data/Sunhh/src/Align/Mugsy/mugsy_x86-64-v1r2.3/mugsyenv.sh 
	"mugsy_exe:s",         # /data/Sunhh/src/Align/Mugsy/mugsy_x86-64-v1r2.3/mugsyWGA
	"mugsy_para:s",        # --distance 2000 --minlength 50 
	                       # --outfile S1_syn --seq S1_all.fitMugsy.fa --aln S1.fitMugsy.maf.xmfa
						   # 1>stdout.S1_syn 2>stderr.S1_syn
################################################################################
HH
	exit 1; 
}

$opts{'cpuN'} //= 1; 
$opts{'seqPerBatch'} //= 1000; 
( defined $opts{'outfile'} and defined $opts{'inmaf'} and defined $opts{'infas'} ) or &usage(); 
$opts{'help'} and &usage(); 
my $HOME = `echo \$HOME`; 
chomp($HOME); $HOME =~ s!\s+$!!; 
$opts{'pl_maf2fasta'}           //= $HOME . '/tools/github/NGS_data_processing/assemble_tools/maf2fasta.pl'; 
$opts{'pl_get_paired_maf'}      //= $HOME . '/tools/github/NGS_data_processing/assemble_tools/get_paired_maf.pl'; 
$opts{'pl_format_maf_forMugsy'} //= $HOME . '/tools/github/NGS_data_processing/assemble_tools/format_maf_forMugsy.pl'; 
$opts{'pl_run_cmd_in_batch'}    //= 'run_cmd_in_batch.pl'; 
$opts{'mugsy_src'} //= '/data/Sunhh/src/Align/Mugsy/mugsy_x86-64-v1r2.3/mugsyenv.sh '; 
$opts{'mugsy_exe'} //= '/data/Sunhh/src/Align/Mugsy/mugsy_x86-64-v1r2.3/mugsyWGA'; 
$opts{'mugsy_para'} //= '--distance 2000 --minlength 50'; 


my $ori_dir = &fileSunhh::_abs_path('./'); 

&tsmsg("[Rec] Beginning in directory: $ori_dir\n"); 

&tsmsg("[Rec] Reading $opts{'inmaf'}\n"); 
my $inFh = &openFH($opts{'inmaf'}, '<'); 
my @all_blks; 
my @all_pairs; 
my %seqId_to_blkIdx; 
while ( my %rec1 = %{ &readMAF($inFh) } ) {
	chomp( $rec1{'a'}[0] ); 
	push(@all_blks, [ $rec1{'a'}[0] ]); 
	my @curr_pairs; 
	for (my $i=0; $i<2; $i++) {
		$rec1{'o'}[$i] =~ m/^s\s/ or &stopErr("Wrong line[$i]: $rec1{o}[$i]\n");
		chomp( $rec1{'o'}[$i] ); 
		my $tr = splitMafSline($rec1{'o'}[$i], 1); 
		push(@curr_pairs, $tr->{'seqId'}); 
		push(@{$all_blks[-1]}, [ $tr->{'seqId'}, $rec1{'o'}[$i] ]); 
		$seqId_to_blkIdx{$tr->{'seqId'}}{ $#all_blks } = 1; 
	}
	push(@all_pairs, [@curr_pairs]); 
}
close($inFh); 

# Sorting. 
&tsmsg("[Rec] Generating groups.\n"); 
my ($grp_aref) = &mathSunhh::divide_group_fromArray(\@all_pairs); 
@$grp_aref = map { [ sort @$_ ] } @$grp_aref; 
my %grp_len = map { $_->[0] => scalar(@$_) } @$grp_aref; 
@$grp_aref = sort { $grp_len{$b->[0]} <=> $grp_len{$a->[0]} || $a cmp $b } @$grp_aref; 
## Recording groups. 
&tsmsg("[Rec] Groups are stored in file [$opts{'outfile'}.grp]\n"); 
my $oGrpFh = &openFH("$opts{'outfile'}.grp", '>'); 
print {$oGrpFh} join("\t", qw/GroupID MemberNum MemberIDs/)."\n"; 
for (my $i=0; $i<@$grp_aref; $i++) {
	print {$oGrpFh} join( "\t", $i, scalar(@{$grp_aref->[$i]}), join(' ', @{$grp_aref->[$i]}) )."\n"; 
}
close ($oGrpFh); 

# Dividing groups and making batch files. 
&tsmsg("[Rec] Separating groups.\n"); 
my $wrk_dir = &fileSunhh::new_tmp_dir(); 
mkdir($wrk_dir); 
chdir($wrk_dir); 
my @new_maf_in; 
my $num_to_reach = $opts{'seqPerBatch'}; 
my $num_have = 0; 
my $tmpFh; 
for (my $i=0; $i<@$grp_aref; $i++) {
	if ($i == 0) {
		my $tmp_file = &fileSunhh::new_tmp_file(); 
		push(@new_maf_in, $tmp_file); 
		$tmpFh = &openFH($tmp_file, '>'); 
		&tell_num(); 
		print {$tmpFh} "##maf version=1\n"; 
	} elsif ( $num_have >= $num_to_reach ) {
		close($tmpFh); 
		while ($num_have >= $num_to_reach) { $num_to_reach += $opts{'seqPerBatch'} } ; 
		my $tmp_file = &fileSunhh::new_tmp_file(); 
		push(@new_maf_in, $tmp_file); 
		$tmpFh = &openFH($tmp_file, '>'); 
		&tell_num(); 
		print {$tmpFh} "##maf version=1\n"; 
	}
	my %has_out; 
	for my $seqId ( @{$grp_aref->[$i]} ) {
		for my $blkIdx ( sort {$a<=>$b} keys %{$seqId_to_blkIdx{$seqId}} ) {
			defined $has_out{$blkIdx} and next; 
			print {$tmpFh} "$all_blks[$blkIdx][0]\n"; 
			print {$tmpFh} "$all_blks[$blkIdx][1][1]\n"; 
			print {$tmpFh} "$all_blks[$blkIdx][2][1]\n"; 
			print {$tmpFh} "\n"; 
			$has_out{$blkIdx} = 1; 
		}
	}
	$num_have += scalar(@{$grp_aref->[$i]}); 
}
close ($tmpFh); 
chdir($ori_dir); 

# Making shell file. 
my $shell_to_run = &fileSunhh::new_tmp_file(); 
&tsmsg("[Rec] Making shell script [$shell_to_run]\n"); 
my $osFh = &openFH($shell_to_run, '>'); 
my @mugsy_raw_out_maf; 
my @mugsy_paired_maf; 
for my $mafName (@new_maf_in) {
	my $shFh = &openFH("$wrk_dir/$mafName.sh", '>'); 
	print {$shFh} "perl $opts{'pl_format_maf_forMugsy'} $wrk_dir/$mafName -specs \",\" -out $wrk_dir/$mafName.fitMugsy\n"; 
	print {$shFh} "perl $opts{'pl_maf2fasta'} $wrk_dir/$mafName.fitMugsy > $wrk_dir/$mafName.xmfa\n"; 
	print {$shFh} "source $opts{'mugsy_src'} ; $opts{'mugsy_exe'} $opts{'mugsy_para'} --outfile $wrk_dir/$mafName.out --seq $opts{infas} --aln $wrk_dir/$mafName.xmfa 1>$wrk_dir/stdout.$mafName 2>$wrk_dir/stderr.$mafName\n"; 
	push(@mugsy_raw_out_maf, "$wrk_dir/$mafName.out.maf"); 
	print {$shFh} "perl $opts{'pl_get_paired_maf'} $wrk_dir/$mafName.out.maf > $wrk_dir/$mafName.paired.maf\n"; 
	push(@mugsy_paired_maf, "$wrk_dir/$mafName.paired.maf"); 
	print {$shFh} "rm $wrk_dir/$mafName.out.maf"; 
	close($shFh); 
	print {$osFh} "$opts{'pl_run_cmd_in_batch'} $wrk_dir/$mafName.sh -cpuN 1 -nprocF $wrk_dir/$mafName.sh.Nproc\n"; 
}
close($osFh); 

# Execute shell script with run_cmd_in_batch.pl 
if ($opts{'printCmd'}) {
	&tsmsg("[ToRUN] $opts{'pl_run_cmd_in_batch'} -cpuN $opts{'cpuN'} $shell_to_run"); 
} else {
	&exeCmd_1cmd("$opts{'pl_run_cmd_in_batch'} -cpuN $opts{'cpuN'} $shell_to_run"); 
	# my $oMergeFh = &openFH("$opts{'outfile'}.raw_Mugsy.maf"); 
	my $oMergeFh = &openFH("$opts{'outfile'}.paired_Mugsy.maf"); 
	my $has_head = 0; 
	for my $mafName (@mugsy_paired_maf) { 
	# for my $mafName (@mugsy_raw_out_maf) {
		my $th = &openFH($mafName,'<'); 
		while (<$th>) {
			if (m/^\s*#/) {
				$has_head == 0 and print {$oMergeFh} $_; 
				$has_head = 1; 
				next; 
			}
			print {$oMergeFh} $_; 
		}
		close($th); 
	}
	close($oMergeFh); 
	# &exeCmd_1cmd( "perl $opts{'pl_get_paired_maf'} $opts{'outfile'}.raw_Mugsy.maf > $opts{'outfile'}.paired.maf" ); 
}

&tsmsg("[Rec] All done for $0.\n"); 

sub tell_num {
	&tsmsg("[Msg] batch_num=", scalar(@new_maf_in)," num_have=$num_have num_to_reach=$num_to_reach\n"); 
}
