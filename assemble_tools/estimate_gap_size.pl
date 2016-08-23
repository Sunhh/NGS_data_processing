#!/usr/bin/perl -w
use strict; 
use warnings; 
use LogInforSunhh; 
use SeqAlnSunhh; 
use mathSunhh; 
use fileSunhh; 


!@ARGV and die "perl $0 in.fa.Nlis pe_aln.srt.bam [max_flank_size min_PE_number]\n"; 

my $fn_nlis = shift; 
my $fn_bam = shift; 
my $max_flank = shift; 
my $min_pe_num = shift; 
$max_flank //= 40000; 
$min_pe_num //= 3; 

my %flag_need_F = %{ &SeqAlnSunhh::mk_flag( 'keep' => '0=1,2=0,3=0,4=0,5=1' , 'drop'=>'' ) }; 
my %flag_need_R = %{ &SeqAlnSunhh::mk_flag( 'keep' => '0=1,2=0,3=0,4=1,5=0' , 'drop'=>'' ) }; 
my %flag_need   = %{ &SeqAlnSunhh::mk_flag( 'keep' => '0=1,2=0,3=0,4=0,5=1;0=1,2=0,3=0,4=1,5=0' , 'drop'=>'' ) }; 

my @nlist = &fileSunhh::load_tabFile( $fn_nlis ); 
#my $fh_nlis = &openFH($fn_nlis, '<'); 
#my %nlist; 
#while (&wantLineC($fh_nlis)) {
#	my @ta=&splitL("\t", $_); 
#	$ta[2] > $ta[3] and @ta[2,3] = @ta[3,2]; 
#	push(@{$nlist{$ta[0]}}, [$ta[2,3]]); 
#}
#close($fh_nlis); 

my %tmp_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>2e3 , 'cntN' => 0 ); 
for my $a1 (@nlist) {
	$tmp_cnt{'cntN'} ++; 
	my ($id, $s, $e) = @{$a1}[0,2,3]; 
	&fileSunhh::log_section( $tmp_cnt{'cntN'} , \%tmp_cnt ) and &tsmsg("[Msg] Processing $tmp_cnt{'cntN'} line. [$id:$s-$e]\n"); 
	if ($id =~ m!^Key$!i) {
		print STDOUT join("\t", @$a1, qw/PE_INSize PE_Num/)."\n"; 
		next; 
	}
	$s > $e and ($s,$e) = ($e,$s); 
	my $s_get = $s-$max_flank; 
	my $e_get = $e+$max_flank; 
	$s_get > 0 or $s_get = 1; 
	my @ins_sizes; 
	open F,'-|',"samtools view $fn_bam $id:${s_get}-${e_get} " or &stopErr("[Err] failed cmd: samtools view $fn_bam $id:${s_get}-${e_get}\n"); 
	my %skip_rd; 
	while (<F>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		defined $skip_rd{$ta[0]} and next; 
		$skip_rd{$ta[0]} = 1; 
		$ta[6] eq '=' or next; 
		defined $flag_need{$ta[1]} or next; 
		my %infor_hash = %{ &SeqAlnSunhh::sam_line2hash( \@ta, [qw/ins_s ins_e/] ) }; 
		$infor_hash{'ins_s'} > $infor_hash{'ins_e'} and next; 
		$infor_hash{'ins_s'} > $e and next; 
		$infor_hash{'ins_e'} < $s and next; 
		push(@ins_sizes, $infor_hash{'ins_e'}-$infor_hash{'ins_s'}+1); 
	}
	close F; 
	my @est_size = ('NA','NA'); 
	if ( @ins_sizes >= $min_pe_num ) {
		my $hr = &mathSunhh::ins_calc( \@ins_sizes, 0 ); 
		$est_size[0] = sprintf("%.1f", $hr->{'interval_mean'}); 
		$est_size[1] = $hr->{'interval_cnt'}; 
	}
	print STDOUT join("\t", @$a1, $est_size[0], $est_size[1])."\n"; 
}


