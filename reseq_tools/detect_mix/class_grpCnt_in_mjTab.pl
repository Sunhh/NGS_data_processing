#!/usr/bin/perl -w
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"ind2grp_list:s", # 
	"colN_start:i",   # 2 
	"cpuN:i", 
	"wind_start:i", 
	"wind_length:i", 
	"wind_step:i", 
	"wind_end:i", 
	"trimTail!", 
); 
$opts{'colN_start'} //= 3; 
$opts{'cpuN'} //= 20; 
$opts{'wind_start'}  //= 1; 
$opts{'wind_length'} //= 1; 
$opts{'wind_step'}   //= $opts{'wind_length'}; 
$opts{'wind_end'}    //= 9999999; 

# Output format : 
# TaxID \\t sPos \\t ePos \\t Grp1ID_Cnt \\t Grp2ID_Cnt \\t Rest_Cnt \\t MISS_Cnt

my $help_txt = <<HH; 
################################################################################
# perl $0 set02.tab.filt.tab_class -ind2grp_list grpMjAL_diff_list.CA_to_CLV > set02.tab.filt.tab_class_byWind
#
# -wind_start       [$opts{'wind_start'}]
# -wind_length      [$opts{'wind_length'}]
# -wind_step        [same to 'wind_length']
# -wind_end         [$opts{'wind_end'}]
#
# -trimTail         [Boolean]
################################################################################
# Example of grpMjAL_diff_list.CA_to_CLV : 
#   CC_CLV     \\t  grp_CLV
#   CC_CLM_CLV \\t  grp_CLV
#   CLM_CLV    \\t  grp_CLV
#   CLV        \\t  grp_CLV
#   CC_CA      \\t  grp_CA
#   CC_CA_CLM  \\t  grp_CA
#   CA_CLM     \\t  grp_CA
#   CA         \\t  grp_CA
################################################################################
# Example of output format : 
#   TaxID \\t chrID \\t windS \\t windE \\t windL \\t grp_CLV_Cnt \\t grp_CA_Cnt \\t Rest_Cnt \\t MISS_Cnt
# Here 'MISS' stands for '.' in './.'; 
# Here 'Rest' includes all tab_classes not defined in -ind2grp_list, including 'REST' class. 
################################################################################
HH


-t and !@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 
defined $opts{'ind2grp_list'} or &LogInforSunhh::usage($help_txt); 

&tsmsg("[Msg] Beigin [$0]\n"); 

my %sample2grp = &load_ind2grpID( $opts{'ind2grp_list'} ); 
my $pm = &LogInforSunhh::get_pm( $opts{'cpuN'} ); 

our @InFp = () ;
if ( !@ARGV )
{
	-t or @InFp = (\*STDIN);
}
else
{
	for (@ARGV) {
		push( @InFp, &openFH($_,'<') );
	}
}


my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
for ( my $i=0; $i<@InFp; $i++ ) {
	my $fh = $InFp[$i]; 
	# header_txt
	my $header_txt = <$fh>; 
	chomp($header_txt); 
	my @ha = &splitL("\t", $header_txt); 
	my $max_cN = $#ha; 
	my @need_keys = ( @{$sample2grp{'grpAr'}}, 'MISS', 'Rest' ); 
	print STDOUT join("\t", 'sampleID', $ha[0], qw/WindS WindE/, @need_keys)."\n"; 
	
	# get sub files
	my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_$i.", 'tmpFile' => "$wrk_dir/base_$i" ); 
	# Process sub files
	&tsmsg("[Msg]   Separating by lines.\n"); 
	for my $sfn (@sub_fn) {
		my $pid = $pm->start and next; 
		open F,'<',"$sfn" or die; 
		my @o_txt; 
		while (<F>) {
			chomp; 
			my @ta = &splitL("\t", $_); 
			for (my $j=$opts{'colN_start'}; $j<@ta; $j++) {
				$o_txt[$j] .= join("\t", @ta[0,1, $j])."\n"; 
				# chr pos CC_CLV/CC_CA_CLV; 
			}
		}
		close F; 
		for (my $j=$opts{'colN_start'}; $j<@ha; $j++) {
			open O,'>',"$sfn.o$j" or die; 
			print O $o_txt[$j]; 
			close O; 
		}
		unlink($sfn); 
		$pm->finish; 
	}
	$pm->wait_all_children; 
	# Merge sub files
	&tsmsg("[Msg]   Separating by columns.\n"); 
	for (my $j=$opts{'colN_start'}; $j<=$max_cN; $j++) {
		my $pid = $pm->start and next; 
		open O,'>',"$wrk_dir/singleC_$j" or die; 
		for my $sfn (@sub_fn) {
			open F,'<',"$sfn.o$j" or die; 
			while (<F>) {
				print O $_; 
			}
			close F; 
			unlink("$sfn.o$j"); 
		}
		close O; 
		$pm->finish; 
	}
	$pm->wait_all_children; 
	# Count windows for each column window. 
	&tsmsg("[Msg]   Computing windows for each column.\n"); 
	for (my $j=$opts{'colN_start'}; $j<=$max_cN; $j++) {
		my $pid = $pm->start and next; 
		open F,'<',"$wrk_dir/singleC_$j" or die; 
		my %chr_wind; 
		while (<F>) {
			chomp; 
			my ($chrID, $chrPos, $tab_cls) = &splitL("\t", $_); 
			$tab_cls =~ m!^([^\s/]+)/([^\s/]+)$! or &stopErr("[Err] Bad tab_class [$tab_cls] : $_\n"); 
			my ($b1, $b2) = ($1, $2); 
			my $g1 = ( $b1 eq '.' ) ? 
			  'MISS' 
			  : ( defined $sample2grp{'grpID'}{$b1} ) ? 
			    $sample2grp{'grpID'}{$b1} 
			    : 'Rest' 
			; 
			my $g2 = ( $b2 eq '.' ) ? 
			  'MISS' 
			  : ( defined $sample2grp{'grpID'}{$b2} ) ? 
			    $sample2grp{'grpID'}{$b2} 
			    : 'Rest' 
			; 
			unless ( defined $chr_wind{$chrID} ) {
				$chr_wind{$chrID} = $ms_obj->setup_windows(
				  'ttl_start' => $opts{'wind_start'}, 
				  'ttl_end'   => $opts{'wind_end'}, 
				  'wind_size' => $opts{'wind_length'}, 
				  'wind_step' => $opts{'wind_step'}, 
				  'minRatio'  => 0
				); 
			}
			my (@wind_i) = @{ $ms_obj->map_windows( 'posi' => $chrPos , 'wind_hash' => $chr_wind{$chrID} ) }; 
			for my $ti (@wind_i) {
				$chr_wind{$chrID}{'cnt'}{$ti}{$g1} ++; 
				$chr_wind{$chrID}{'cnt'}{$ti}{$g2} ++; 
			}
		}
		my %endIdx; 
		if ( $opts{'trimTail'} ) {
			for my $chrID ( sort keys %chr_wind ) {
				my @idx = @{ $chr_wind{$chrID}{'info'}{'windSloci'} }; 
				my $last_i = $#idx; 
				for ( ; $last_i >= 0; $last_i -- ) {
					defined $chr_wind{$chrID}{'cnt'}{$idx[$last_i]} and last; 
				}
				$endIdx{$chrID} = ( $last_i == -1 ) ? 'stop' : $idx[$last_i]; 
			}
		}
		open O,'>',"$wrk_dir/singleCWind_$j" or die; 
		for my $chrID ( sort keys %chr_wind ) {
			for my $ti ( @{$chr_wind{$chrID}{'info'}{'windSloci'}} ) {
				for my $k1 ( @need_keys ) {
					$chr_wind{$chrID}{'cnt'}{$ti}{$k1} //= 0; 
				}
				$opts{'trimTail'} and ( $endIdx{$chrID} eq 'stop' or $endIdx{$chrID} < $ti ) and last; 
				print O join("\t", $ha[$j], $chrID, @{$chr_wind{$chrID}{'loci'}{$ti}}[0,1], @{$chr_wind{$chrID}{'cnt'}{$ti}}{@need_keys})."\n"; 
			}
		}
		close O; 
		close F; 
		$pm->finish; 
	}
	$pm->wait_all_children; 
	# Combine each file "$wrd_dir/singleCWind_$j"; 
	&tsmsg("[Msg]   Combine all singleCWind_ files\n"); 
	for (my $j=$opts{'colN_start'}; $j<=$max_cN; $j++) {
		open F,'<',"$wrk_dir/singleCWind_$j" or die; 
		while (<F>) {
			print STDOUT $_; 
		}
		close F; 
	}
	# Delete temporary files 
	&fileSunhh::_rmtree($wrk_dir); 
}

&tsmsg("[Msg] Finish [$0]\n"); 



sub load_ind2grpID {
### Input format : 
# sample1  \\t   grp3
# sample2  \\t   grp3
# sample3  \\t   grp1
### Output format : 
# $h{'order'}{'sample1'} = 1
# $h{'order'}{'sample2'} = 1
# $h{'order'}{'sample3'} = 2
# $h{'grpID'}{'sample1'} = 'grp3'
# $h{'grpID'}{'sample2'} = 'grp3'
# $h{'grpID'}{'sample3'} = 'grp1'
# $h{'grpAr'} = [ 'grp3', 'grp1' ]
	my %h;
	my %reset_num;
	my $reset_num_next = 1;
	my $fh = &openFH( $_[0], '<' ); 
	while (<$fh>) {
		chomp;
		my @ta = split(/\t/, $_);
		unless ( defined $reset_num{$ta[1]} ) {
			$reset_num{$ta[1]} = $reset_num_next;
			$reset_num_next ++;
			push(@{$h{'grpAr'}}, $ta[1]); 
		}
		$h{'order'}{$ta[0]} = $reset_num{$ta[1]}; 
		$h{'grpID'}{$ta[0]} = $ta[1]; 
	}
	close($fh);
	return %h;
}# load_ind2grpID ()


