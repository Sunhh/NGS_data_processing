#!/usr/bin/perl
# 2020-03-31 Adapt to repeatmodeler 2.0; 
use strict; 
use warnings; 
use LogInforSunhh; 
use Parallel::ForkManager;
use File::Path qw(make_path remove_tree);
use File::Copy; 
use Cwd 'abs_path'; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 
use fileSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"path_repClass:s", "para_RC:s", 
	"inFa:s", "tagUnknown:s", 
	"noChopID!", 
	"cpuN:i", 
	"out:s", "outTbl:s", 
	"repeat_times:i", 
	"v2_modeler!", 
); 

$opts{'repeat_times'} //= 1; 
$opts{'path_repClass'} = $opts{'path_repClass'} // '/share/app/Annotation/repeatmodeler/RepeatModeler/RepeatClassifier'; 
$opts{'para_RC'} = $opts{'para_RC'} // '-engine ncbi'; 
$opts{'tagUnknown'} = $opts{'tagUnknown'} // 'LTR'; 
$opts{'cpuN'} = $opts{'cpuN'} // 1; 
$opts{'pl_dealFa'} = $opts{'pl_dealFa'} // `echo \$HOME/tools/github/NGS_data_processing/deal_fasta.pl`; 
chomp( $opts{'pl_dealFa'} ); 

$opts{'help'} and &usage(); 
defined $opts{'inFa'} and -e $opts{'inFa'} or &usage(); 
$opts{'out'} = $opts{'out'} // "$opts{'inFa'}.classified"; 
$opts{'outTbl'} = $opts{'outTbl'} // "$opts{'out'}.tbl"; 

&tsmsg("[Rec] Begin $0\n"); 

my $wrk_dir = &fileSunhh::new_tmp_dir('create' => 1); 
my $cur_dir = abs_path( File::Path::getcwd() ); 
my $tmp_file = "tmp_toClass.fa"; 
my %ori_seq = %{$fs_obj->save_seq_to_hash( 'faFile' => $opts{'inFa'}, 'has_head' => 1 )}; 
my %ori_cls; 
my %got_cls; 
open OT,'>',"$wrk_dir/$tmp_file" or die; 
for my $k1 (sort { $ori_seq{$a}{'Order'} <=> $ori_seq{$b}{'Order'} } keys %ori_seq) {
	my %t1 = %{ &parse_input_header( $ori_seq{$k1}{'head'} ) }; 
	$ori_cls{ $t1{'key_shrt'} } = $t1{'class'}; 
	my $oheader = $t1{'header_shrt'}; 
	if ($opts{'v2_modeler'}) {
		$oheader .= "#Unknown"; 
	}

	chomp( $ori_seq{$k1}{'seq'} ); 	
	print OT ">$oheader\n$ori_seq{$k1}{'seq'}\n"; 
}
close OT; 
chdir($wrk_dir); 
if ( $opts{'cpuN'} <= 1 ) {
	for (my $i=0; $i<$opts{'repeat_times'}; $i++) {
		for my $fn ('tmpConsensi.fa.cat.all', "$tmp_file.classified") {
			unlink($fn); 
		}
		&exeCmd("$opts{'path_repClass'} -consensi $tmp_file $opts{'para_RC'}"); 
		&parse_cls_file( "$tmp_file.classified", "$tmp_file.cls_tbl", '>>' ); 
		%got_cls = %{ &sum_cls_tbl( "$tmp_file.cls_tbl" ) }; 
	}
} else {
	
	&exeCmd("perl $opts{'pl_dealFa'} $tmp_file -cut 1 -cut_prefix cut1 1>cut1.std 2>cut1.err"); 
	mkdir("run_cut1"); 
	opendir DD,"cut1_cutted/" or &stopErr("[Err] $!\n");
	my @inFas; 
	for my $fn (readdir(DD)) {
		$fn =~ m/^\./ and next; 
		$fn =~ m/^cut1_\d+.fasta$/ or next; 
		mkdir("run_cut1/dir_${fn}/", 0755); 
		File::Copy::copy("cut1_cutted/$fn", "run_cut1/dir_${fn}/$fn") or &stopErr("[Err] Failed to copy file $fn\n"); 
		push(@inFas, $fn); 
	}
	closedir(DD); 
	@inFas = sort { my ($n1) = ($a =~ m/_(\d+)\.fasta$/); my ($n2) = ($b =~ m/_(\d+)\.fasta$/); $n1<=>$n2; } @inFas; 
	
	my $cwd = abs_path( File::Path::getcwd() ); 
	
	my $MAX_PROCESSES = $opts{'cpuN'} ; # Sometimes $parm{'cpuN'} - 1 may be better.
	my $pm = new Parallel::ForkManager($MAX_PROCESSES); 
	for my $fn (@inFas) {
		my $pid = $pm->start and next; 
		chdir("run_cut1/dir_${fn}/"); 
		for (my $i=0; $i<$opts{'repeat_times'}; $i++) {
			for my $tfn ('tmpConsensi.fa.cat.all', "${fn}.classified") {
				unlink($tfn); 
			}
			&exeCmd("$opts{'path_repClass'} -consensi $fn $opts{'para_RC'}"); 
			&parse_cls_file( "${fn}.classified", "${fn}.cls_tbl", '>>' ); 
		}
		chdir("../../"); 
		$pm->finish; 
	}
	$pm->wait_all_children; 
	
	chdir($cwd); 
	-e "${tmp_file}.cls_tbl" and unlink("${tmp_file}.cls_tbl"); 
	for my $fn (@inFas) {
		-e "run_cut1/dir_${fn}/${fn}.cls_tbl" or do { &tsmsg("[Err] File run_cut1/dir_${fn}/${fn}.cls_tbl not found.\n"); next; }; 
		&exeCmd("cat run_cut1/dir_${fn}/${fn}.cls_tbl >> ${tmp_file}.cls_tbl"); 
	}
	%got_cls = %{ &sum_cls_tbl( "$tmp_file.cls_tbl" ) }; 
}
chdir($cur_dir); 

open BB, '>',"$opts{'out'}"    or &stopErr("[Err] $!\n"); 
open BBT,'>',"$opts{'outTbl'}" or &stopErr("[Err] $!\n"); 
print BBT join("\t",qw/Raw_ID Raw_Class New_Class Got_CLass/)."\n"; 
for my $k1 (sort { $ori_seq{$a}{'Order'} <=> $ori_seq{$b}{'Order'} } keys %ori_seq) {
	my %t1 = %{ &parse_input_header( $ori_seq{$k1}{'head'} ) }; 
	my $key_shrt = $t1{'key_shrt'}; 
	defined $got_cls{ $key_shrt } or &stopErr("[Err] Failed to get class for [$k1]\n"); 
	my @cls_arr = map { [ $_, $got_cls{$key_shrt}{'Class'}{$_} ] } sort { $got_cls{$key_shrt}{'Class'}{$b} <=> $got_cls{$key_shrt}{'Class'}{$a} || $got_cls{$key_shrt}{'ClassIdx'}{$a} <=> $got_cls{$key_shrt}{'ClassIdx'}{$b} } keys %{ $got_cls{$key_shrt}{'Class'} }; 
	my $cls_txt = join(';;', map { "$_->[0]:$_->[1]" } @cls_arr); 
	my $cls_use = $cls_arr[0][0]; 
	if ($cls_use =~ m/^Unknown$/i) {
		if ( defined $ori_cls{$key_shrt} and $ori_cls{$key_shrt} !~ m/^Unknown$/i and $ori_cls{$key_shrt} ne '' ) {
			$cls_use = $ori_cls{$key_shrt}; 
		} else {
			$cls_use = $opts{'tagUnknown'}; 
		}
	}
	my $h1 = "$t1{'key_shrt'}#$cls_use$t1{'rest_header'}"; 
	print BB ">$h1\n$ori_seq{$k1}{'seq'}\n"; 
	print BBT join("\t", $k1, $t1{'class'}, $cls_use, $cls_txt)."\n"; 
}
close BBT; 
close BB; 

&fileSunhh::_rmtree($wrk_dir); 

&tsmsg("[Rec] Finish $0\n"); 

sub usage {
	print STDOUT <<HH; 
################################################################################
# perl $0 -inFa toBeClassified.fas -tagUnknown LTR
#
# -help
# -path_repClass [RepeatClassifier]
# -v2_modeler    [Boolean] Give this if using RepeatModeler V2.0; 
# -para_RC       [-engine ncbi]
#
# -out           [inFa.classified]
# -outTbl        [inFa.classified.tbl]
#
# -noChopID      [boolean] Do not chop ID if given. 
#
# -cpuN          [1] Parallelly run if larger than 1. 
# -pl_dealFa     [deal_fasta.pl] Must be given if cpuN > 1
#
# -repeat_times  [$opts{'repeat_times'}]
#
################################################################################
HH
	exit 1; 
}
sub parse_cls_file {
	my ($fin, $fout, $ow) = @_; 
	$fout //= "${fin}.cls_tbl"; 
	$ow   //= '>'; 
	my %ss = %{ $fs_obj->save_seq_to_hash('faFile'=>$fin, 'has_head'=>1) }; 
	my $ofh = &openFH($fout, $ow); 
	for my $k0 (sort {$ss{$a}{'Order'} <=> $ss{$b}{'Order'}} keys %ss) {
		$k0 =~ m!^(\S+)\#([^#\s]+)$! or &stopErr("[Err] Failed to parse ID [$k0]\n"); 
		print {$ofh} join("\t", $1, $2)."\n"; 
	}
	close ($ofh); 
	return; 
}# parse_cls_file() 
sub sum_cls_tbl {
	my ($fin) = @_; 
	my %back; 
	open F,'<',"$fin" or &stopErr("[Err] Failed [$fin] $!\n"); 
	while (<F>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		$back{$ta[0]}{'Order'} //= $.; 
		$back{$ta[0]}{'Class'}{$ta[1]} ++; 
		$back{$ta[0]}{'ClassIdx'}{$ta[1]} //= $.; 
	}
	close F; 
	return(\%back); 
}# sum_cls_tbl() 

sub parse_input_header {
	my ($t_header) = @_; 
	chomp($t_header); 
	my %back; 
	$t_header =~ s!^(\S+)!! or &stopErr("[Err] Bad header [$t_header]\n"); 
	$back{'rest_header'} = $t_header; 
	my $t_key = $1; 
	$back{'key_long'} = $t_key; 
	if ( $t_key =~ m!^(\S+)\#([^\s#]+)$! ) {
		$back{'key_shrt'} = $1; 
		$back{'class'}    = $2; 
	} else {
		$back{'key_shrt'} = $1; 
		$back{'class'}    = ''; 
	}
	unless ( $opts{'noChopID'} ) {
		$t_header = $back{'key_shrt'}; 
	}
	$back{'header_shrt'} = $t_header; 
	return(\%back); 
} # parse_input_header () 

