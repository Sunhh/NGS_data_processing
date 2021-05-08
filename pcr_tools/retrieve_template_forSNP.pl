#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"l_flank_size:i", # 300
	"r_flank_size:i", # 300
	"ref_fa:s",       # /Data/Sunhh/melon/db/melon_v351edit.chr.fa
	"out_prefix:s",   # out
	"min_flank_size:i", # 50 
); 

$opts{'l_flank_size'} //= 300; 
$opts{'r_flank_size'} //= 300; 
$opts{'ref_fa'}       //= '/Data/Sunhh/melon/db/melon_v351edit.chr.fa'; 
$opts{'min_flank_size'} //= 50; 

my $htxt = <<HHH; 
################################################################################
# perl $0 site.list -out_prefix opre 
# 
# -ref_fa         [fasta file path] $opts{'ref_fa'}
# -min_flank_size [Num] $opts{'min_flank_size'}
# 
# -help
################################################################################
#  Example of site.list: 
#  Site_ID             chr_ID  chr_Position    Base_Ref           Base_indv1             Base_indv2
#  Site.HS200616.01    chr1    3605732         T                  TAGACTTTACTAAACGATC    T
#  Site.HS200616.02    chr3    3103223         T                  TCTAAAAAGTATTCATGTA    T
#  Site.HS200616.03    chr3    4583857         GTTTAAATGTTTTAGCA  G                      GTTTAAATGTTTTAGCA
################################################################################
HHH
$opts{'help'} and &LogInforSunhh::usage($htxt); 
-t and !@ARGV and &LogInforSunhh::usage($htxt); 

my %seqs = %{ $fs_obj->save_seq_to_hash('faFile'=> $opts{'ref_fa'}) }; 
for (keys %seqs) { $seqs{$_}{'seq'} =~ s!\s!!g; $seqs{$_}{'len'} = length($seqs{$_}{'seq'});  }

my $of1 = "$opts{'out_prefix'}.temp1.fas"; &fileSunhh::write2file($of1, '', '>'); 
my $of2 = "$opts{'out_prefix'}.tempX.tab"; &fileSunhh::write2file($of2, '', '>'); 
&fileSunhh::write2file($of2, join("\t", qw/Site_ID targetSE.temp1 len.flankLR loc.temp0 seq.temp1 seq.temp2 seq.temp0/)."\n", '>>'); 

while (<>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	$ta[0] eq 'Site_ID' and $ta[1] eq 'chr_ID' and next; 
	defined $seqs{$ta[1]} or do { &tsmsg("[Wrn] No chrID [$ta[0]] found.\n"); next; }; 
	my $ins_0 = length($ta[3]); 
	my $ins_1 = length($ta[4]); 
	my $ins_2 = length($ta[5]); 
	my ($l_s, $l_e) = ( $ta[2]-$opts{'l_flank_size'}, $ta[2]-1 ); 
	my ($r_s, $r_e) = ( $ta[2]+$ins_0-1+1, $ta[2]+$ins_0-1+$opts{'r_flank_size'} ); 
	$l_s < 1 and $l_s = 1; 
	$r_e > $seqs{$ta[1]}{'len'} and $r_e = $seqs{$ta[1]}{'len'}; 
	my $l_seq = substr($seqs{$ta[1]}{'seq'}, $l_s-1, $l_e-$l_s+1); 
	# $l_seq =~ m![nNxX]$! and do { &tsmsg("[Wrn] Not enough left flank for line: $_\n"); next; }; 
	$l_seq =~ s!^(.*[nNxX]+)([^nNxX]*)$!$2! and $l_s += length($1); 
	$l_e-$l_s+1 < $opts{'min_flank_size'} and do { &tsmsg("[Wrn] Not enough left flank for line: $_\n"); next; }; 
	my $r_seq = substr($seqs{$ta[1]}{'seq'}, $r_s-1, $r_e-$r_s+1); 
	# $r_seq =~ m!^[nNxX]! and do { &tsmsg("[Wrn] Not enough right flank for line: $_\n"); next; }; 
	$r_seq =~ s!^([^nNxX]+)([nNxX].*)$!$1! and $r_e -= length($2); 
	$r_e-$r_s+1 < $opts{'min_flank_size'} and do { &tsmsg("[Wrn] Not enough right flank for line: $_\n"); next; }; 
	my $temp0_seq = join('', $l_seq, $ta[3], $r_seq); 
	my $temp1_seq = join('', $l_seq, $ta[4], $r_seq); 
	my $temp2_seq = join('', $l_seq, $ta[5], $r_seq); 
	my $siteS_in_temp1 = $l_e-$l_s+1+1; 
	my $siteE_in_temp1 = $siteS_in_temp1 + $ins_1 - 1; 
	&fileSunhh::write2file($of1, ">$ta[0].temp1\n$temp1_seq\n", '>>'); 
	&fileSunhh::write2file($of2, 
		join("\t", 
			$ta[0], 
			"${siteS_in_temp1}-${siteE_in_temp1}", 
			join(":", $l_e-$l_s+1, $r_e-$r_s+1), 
			join(":", $ta[1], "$l_s-$l_e", join("-",$l_e+1,$r_s-1), "$r_s-$r_e"), 
			$temp1_seq, 
			$temp2_seq, 
			$temp0_seq
		)."\n", 
		'>>'
	); 
}

