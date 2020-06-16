#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"in_tempTab:s", # site.tempX.tab, output of retrieve_template_forSNP.pl; 
	"pl_run_primer3_general:s", # /home/Sunhh/tools/github/NGS_data_processing/pcr_tools/run_primer3_general.pl
	"out_prefix:s", # opref
	"sconf_primer_size:s",   # [opt:min:max] from run_primer3_general.pl 
	"sconf_product_size:s@", # from run_primer3_general.pl 
	"sconf_return_num:i",    # from run_primer3_general.pl
); 

$opts{'out_prefix'}             //= 'opref'; 

$opts{'pl_run_primer3_general'} //= '/home/Sunhh/tools/github/NGS_data_processing/pcr_tools/run_primer3_general.pl'; 
$opts{'sconf_primer_size'}      //= '20:18:27'; # [opt:min:max] from run_primer3_general.pl 
$opts{'sconf_product_size'}     //= [ '151-200', '201-250', '101-150', '251-350' ]; # from run_primer3_general.pl 
my $gg_product_size  = join(" ", map { "-sconf_product_size $_" } @{$opts{'sconf_product_size'}}); 
$opts{'sconf_return_num'}       //= 10; 


my $htxt = <<HHH; 
perl $0 -in_tempTab site.tempX.tab -out_prefix opref
HHH


$opts{'help'} and &LogInforSunhh::usage($htxt); 
defined $opts{'in_tempTab'} or &LogInforSunhh::usage($htxt); 


open F1,'<',"$opts{'in_tempTab'}" or die; 
my $wdir = &fileSunhh::new_tmp_dir('create' => 1); 
open O1,'>',"$opts{'out_prefix'}.primer.tab" or die; 
my %o_hash; 
while (<F1>) {
	chomp; 
	my @ta=split(/\t/, $_); # 
	$ta[0] eq 'Site_ID' and next; 
	# 0       Site_ID Site.HS200616.01
	# 1       targetSE.temp1  301-319
	# 2       len.flankLR     300:300
	# 3       loc.temp0       chr1:3605432-3605731:3605732-3605732:3605733-3606032
	# 4       seq.temp1       CTTTCAT
	# 5       seq.temp2       CTTTCAT
	# 6       seq.temp0       CTTTCAT
	
	$ta[3] =~ m!^(\S+):(\d+)\-(\d+):(\d+)\-(\d+):(\d+)\-(\d+)$! or die "Bad input: [$ta[3]]: $_\n"; 
	my ($chr_ID, $chr_ls, $chr_le, $chr_is, $chr_ie, $chr_rs, $chr_re) = ($1,$2,$3,$4,$5,$6,$7); 

	# Prepare -in_seq; 
	&fileSunhh::write2file("$wdir/in_seq.fas", ">$ta[0].temp1\n$ta[4]\n", '>'); 

	# Prepare -in_loc; 
	my $locID = $ta[0]; 
	my $seqID = "$ta[0].temp1"; 
	$ta[1] =~ m!^(\d+)\-(\d+)$! or do { &tsmsg("[Wrn] Failed to parse targetSE.temp1 [$ta[1]]\n"); next; }; 
	my ($target_start, $target_end) = ($1, $2); 
	my $temp_start = 1; 
	my $temp_end   = length($ta[4]); 
	&fileSunhh::write2file("$wdir/in_loc.tab", 
		join("\t", 
			$locID, 
			$seqID, 
			$target_start, 
			$target_end, 
			$temp_start, 
			$temp_end
		)."\n", 
		'>'
	); 

	# Prepare -in_p3conf
	### No need currently. 
	&run_cmd("perl $opts{'pl_run_primer3_general'} $gg_product_size -sconf_primer_size $opts{'sconf_primer_size'} -sconf_return_num $opts{'sconf_return_num'} -in_seq $wdir/in_seq.fas -in_loc $wdir/in_loc.tab > $wdir/out.tab"); 
	open F2,'<',"$wdir/out.tab" or die; 
	while (my $l_f2 = <F2>) {
		chomp($l_f2); 
		my @tb = split(/\t/, $l_f2); # 
		if ( $tb[9] ne 'l_5pPos' ) {
			$tb[8] = $chr_ID; 
			$tb[9] = $chr_ls + $tb[9] -1; 
			$tb[10] = $chr_ls + $tb[10] -1; 
			$tb[11] = $chr_re - ($temp_end-$tb[11]+1) + 1; 
			$tb[12] = $chr_re - ($temp_end-$tb[12]+1) + 1; 
		}
		my $o_line = join("\t", @tb); 
		defined $o_hash{$o_line} or print O1 "$o_line\n"; 
		$o_hash{$o_line} = 1; 
	}
	close F2; 
}
close F1; 
close O1; 
&fileSunhh::_rmtree($wdir); 

sub run_cmd {
	&exeCmd_1cmd("$_[0]") and &stopErr("[Err] Failed at CMD: $_[0]\n"); 
}



