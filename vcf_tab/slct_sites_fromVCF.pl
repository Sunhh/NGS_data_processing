#!/usr/bin/perl -w
use strict; 
use warnings; 
use Getopt::Long; 
use LogInforSunhh; 
use fileSunhh; 
use fastaSunhh; 
my $fas_obj = fastaSunhh->new(); 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"in_vcf:s",    # Ari_126groups_genotyping.recode.vcf 
	"in_keyLen:s", # db/ITAG2.3_genomic.fa.key_len . Format: key \\t len
	"in_region:s", # Ari_126groups_need_region_tab . Format: lineID \\t grpID \\t  grpFamID \\t targetRegion \\n
	               #                                         3882   \\t 1     \\t  3882     \\t SL2.40ch2:46M-end;;SL2.40ch5:1M-3M;;SL2.40ch6:45M-end;;SL2.40chr7:start-end
	"ref_fas:s",   # /data/qiyue/database/db_bowtie2/ITAG2.3_genomic.fa 
	"out_vcf:s",   # 

	"jar_gatk:s",  # /data/qiyue/tools/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar
	"jar_exe:s",   # java 
); 

$opts{'jar_gatk'} //= '/data/qiyue/tools/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar'; 
$opts{'jar_exe'}  //= 'java'; 

my $help_txt = <<HH; 
####################################################################################################
# perl $0    -in_vcf  Ari_126groups_genotyping.recode.vcf   -in_keyLen db/ITAG2.3_genomic.fa.key_len  -in_region Ari_126groups_need_region_tab  -ref_fas /data/qiyue/database/db_bowtie2/ITAG2.3_genomic.fa -out_vcf slct.vcf
#
# -help 
#
# -in_vcf       [filename] GATK_output.vcf
# -in_keyLen    [filename] Format: key \\t len 
# -in_region    [filename] Format: lineID \\t grpID \\t  grpFamID \\t targetRegion
#                                  3882   \\t 1     \\t  3882     \\t SL2.40ch2:46M-end;;SL2.40ch5:1M-3M;;SL2.40ch6:45M-end;;SL2.40chr7:start-end
#
# -out_vcf      [filename] output_slct.vcf
# 
# -ref_fas      [filename] Should be indexed with .dict and .fai files. 
#
# -jar_exe      [$opts{'jar_exe'}]
# -jar_gatk     [$opts{'jar_gatk'}]
#
####################################################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
for my $t1 (qw/out_vcf in_vcf in_region/) {
	defined $opts{$t1} or &LogInforSunhh::usage($help_txt); 
}

open OV,'>',"$opts{'out_vcf'}" or &stopErr("[Err] Failed to write file [$opts{'out_vcf'}]\n"); 
close OV; 

my %chrLen = %{ &load_keyLen($opts{'in_keyLen'}) }; 
my @tgtLoc = @{ &load_region_tab( $opts{'in_region'}, \%chrLen ) }; 

my $tmp_dir = &fileSunhh::new_tmp_dir(); 
mkdir($tmp_dir) or &stopErr("[Err] Failed to create dir [$tmp_dir/]\n"); 

for my $tgt (@tgtLoc) {
	my ($grpFamID, $locAR) = @$tgt; 
	&tsmsg("[Msg]   Processing grpFamID [$grpFamID]\n"); 
	if ( @$locAR == 0 ) {
		&exeCmd_1cmd("$opts{'jar_exe'} -jar $opts{'jar_gatk'} -T SelectVariants -R $opts{'ref_fas'} -V $opts{'in_vcf'} -o $tmp_dir/curr.vcf -sn $grpFamID") and &stopErr("[Err] CMD failed.\n"); 
	} else {
		open O,'>',"$tmp_dir/curr.list" or &stopErr("[Err] Failed to open [$tmp_dir/curr.list]\n"); 
		for my $loc (@$locAR) {
			print O "$loc\n"; 
		}
		close O; 
		&exeCmd_1cmd("$opts{'jar_exe'} -jar $opts{'jar_gatk'} -T SelectVariants -R $opts{'ref_fas'} -V $opts{'in_vcf'} -o $tmp_dir/curr.vcf -sn $grpFamID -L $tmp_dir/curr.list") and &stopErr("[Err] CMD failed.\n"); 
	}

	if ( -e "$tmp_dir/prev.vcf" ) {
		&exeCmd_1cmd("$opts{'jar_exe'} -jar $opts{'jar_gatk'} -T CombineVariants -R $opts{'ref_fas'} --variant:foo $tmp_dir/prev.vcf --variant:bar $tmp_dir/curr.vcf -o $tmp_dir/comb.vcf -genotypeMergeOptions PRIORITIZE -priority foo,bar") and &stopErr("[Err] CMD failed.\n"); 
		&fileSunhh::_move( "$tmp_dir/comb.vcf", "$tmp_dir/prev.vcf" ); 
	} else {
		&fileSunhh::_move( "$tmp_dir/curr.vcf", "$tmp_dir/prev.vcf" ); 
	}
}

&fileSunhh::_move( "$tmp_dir/prev.vcf", "$opts{'out_vcf'}" ); 

&tsmsg("[Rec] All done [$0].\n"); 

sub load_region_tab {
# lineID  grpID   grpFamID    targetRegion
# 3882    1       3882        SL2.40ch2:46M-end;;SL2.40ch5:1M-3M;;SL2.40ch6:45M-end
# 3886    2       3886        SL2.40ch7:62M-end;;SL2.40ch11:1M-3M;;SL2.40ch11:53M-end
# 3892    3       3892        SL2.40ch11:53M-end
# 1234    4       1234-i      all
	my ($fn, $kl_href) = @_; 
	my $fh = &openFH($fn, '<'); 
	my @back; 
	my %h; 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		my ($grpFamID, $target) = @ta[2,3]; 
		$grpFamID =~ m/^grpFamID$/i and next; 
		defined $h{$grpFamID} and die "repeated grpFamID [$grpFamID]\n"; 
		$h{$grpFamID} = 1; 
		my @tb; 
		if ( $target eq '' or $target =~ m!^\s*all\s*$!i) {
			push(@back, [ $grpFamID, [] ]); 
			next; 
		}
		for my $loc ( &splitL(";;", $target) ) {
			$loc =~ m!^(\S+):(\S+)\-(\S+)$!i or die "loc=|$loc|\n"; 
			my ($cid, $cs, $ce) = ( $1, $2, $3 ); 
# $cs eq 'end' and do { defined $chrLen{$cid} or die "$loc\n"; }; 
# $ce eq 'end' and do { defined $chrLen{$cid} or die "$loc\n"; }; 
			$cs = &get_pos( $cs, $chrLen{$cid} ); 
			$cs == 0 and $cs = 1; 
			$ce = &get_pos( $ce, $chrLen{$cid} ); 
			if ( defined $chrLen{$cid} and $chrLen{$cid} > 0 ) { $ce > $chrLen{$cid} and $ce = $chrLen{$cid}; }
			push(@tb, "$cid:$cs-$ce"); 
		}
		push(@back, [ $grpFamID, [@tb] ]); 
	}
	close($fh); 
	return(\@back); 
}# load_region_tab() 

sub get_pos {
	my ($p, $max) = @_; 
	if ( $p =~ m!^\s*(\d+)\s*$!i ) {
		return $1; 
	} elsif ( $p =~ m!^\s*([\d.]+)\s*K\s*$!i ) {
		return $1*1e3; 
	} elsif ( $p =~ m!^\s*([\d.]+)\s*M\s*$!i ) {
		return $1*1e6; 
	} elsif ( $p =~ m!^\s*([\d.]+)\s*G\s*$!i ) {
		return $1*1e9; 
	} elsif ( $p =~ m!^\s*start\s*$!i ) {
		return 1; 
	} elsif ( $p =~ m!^\s*end\s*$!i ) {
		defined $max or &stopErr("[Err] No max defined.\n"); 
		return $max; 
	} else {
		&stopErr("[Err] Faield to parse position [$p]\n"); 
	}
	return; 
}# get_pos() 

sub load_keyLen {
	my $fn = shift; 
	my %back; 
	unless ( defined $fn ) {
		my %h = %{ $fas_obj->save_seq_to_hash( 'faFile' => $opts{'ref_fas'} ) }; 
		for my $k (keys %h) {
			$h{$k}{'seq'} =~ s!\s!!g; 
			$back{$k} = length( $h{$k}{'seq'} ); 
		}
		return(\%back); 
	}
	my $fh = &openFH($fn, '<'); 
	while (&wantLineC($fh)) {
		my @ta=&splitL("\t", $_); 
		$ta[0] =~ m/^key$/i and next; 
		$back{$ta[0]} = $ta[1]; 
	}
	close ($fh); 
	return(\%back); 
}# load_keyLen()

