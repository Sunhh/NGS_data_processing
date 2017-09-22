#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use SNP_tbl; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"startColN:i", # 2 
	"noHeader!", 
	"cpuN:i", # 20 
	"input_is_tab!", 
); 
$opts{'startColN'} //= 2; 
$opts{'cpuN'} //= 20; 
$opts{'cpuN'} < 1 and $opts{'cpuN'} = 1; 

my $pm = &LogInforSunhh::get_pm( $opts{'cpuN'} ); 
my %d2b_list = &SNP_tbl::get_diploid_d2b(); 
my %used; 
$used{'bad_geno'} = {}; 

my $help_txt = <<HH; 

Count N (missing), heterozygous, and homozygous sites per individual. 
The output format is : qw/IndvID N_Num Typed_Num Het_Num Hom_Num Het_Ratio Hom_Ratio N_Ratio/

perl $0 in.snp > in.snp.cntNHH

 -help               

 -cpuN               [$opts{'cpuN'}]
 -startColN          [$opts{'startColN'}]
 -noHeader           [Bool]

 -input_is_tab       [Bool] Should be given if input is .vcf.tab format. 

Genotype column start from colN=$opts{'startColN'}
Do not parse the first line. 

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my $fh = \*STDOUT; 
@ARGV > 0 and $fh = &openFH( $ARGV[0], '<' ); 

my @ha; 
unless ($opts{'noHeader'}) {
	my $head = <$fh>; 
	chomp($head); 
	@ha=split(/\t/, $head); 
}

# separate files
my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 );
my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" );
close($fh); 
if ( $opts{'noHeader'} and @ha == 0 ) {
	open F,'<',"$sub_fn[0]" or die; 
	$_ = <F>; 
	chomp; 
	my @ta = &splitL("\t", $_); 
	for (my $i=0; $i<@ta; $i++) {
		$ha[$i] = "Col_$i"; 
	}
	close F; 
}

# Process sub-files
for my $sfn ( @sub_fn ) {
        my $pid = $pm->start and next;
        open F,'<',"$sfn" or &stopErr("[Err] Failed to open subfile [$sfn]\n");
        open O,'>',"$sfn.o" or &stopErr("[Err] Failed to write subfile [$sfn.o]\n");
	my @sub_cnt; 
        while (<F>) {
                chomp;
                my @ta = &splitL("\t", $_);
		my @tb = @ta[ $opts{'startColN'} .. $#ta ]; 
		$opts{'input_is_tab'} or &SNP_tbl::aref_cols2tab( \@tb,  $used{'bad_geno'}, \%d2b_list, 0); 
		for (my $i=0; $i<@tb; $i++) {
			my $j=$i+$opts{'startColN'}; 
			my $type = &type_tabGeno( $tb[$i] ); 
			$sub_cnt[$j]{$type} ++; 
		}
        }
	print O join("\t", qw/IndvID N_Num Typed_Num Het_Num Hom_Num Het_Ratio Hom_Ratio N_Ratio/)."\n"; 
	for (my $i=$opts{'startColN'}; $i<@ha; $i++) { 
		for my $tk (qw/miss homo hete/) {
			$sub_cnt[$i]{$tk} //= 0; 
		}
		my $tot        = $sub_cnt[$i]{'miss'} + $sub_cnt[$i]{'hete'} + $sub_cnt[$i]{'homo'}; 
		my $tot_typed  = $sub_cnt[$i]{'hete'} + $sub_cnt[$i]{'homo'}; 
		my $rat_hete   = ($tot_typed > 0) ? ( $sub_cnt[$i]{'hete'}/$tot_typed*100 ) : -1; 
		my $rat_homo   = ($tot_typed > 0) ? ( $sub_cnt[$i]{'homo'}/$tot_typed*100 ) : -1; 
		my $rat_miss   = ($tot       > 0) ? ( $sub_cnt[$i]{'miss'}/$tot      *100 ) : -1; 
		print O join("\t", $ha[$i], $sub_cnt[$i]{'miss'}, $tot_typed, @{$sub_cnt[$i]}{qw/hete homo/}, $rat_hete, $rat_homo, $rat_miss)."\n"; 
	}
        close O;
        close F;
        $pm->finish;
}

$pm->wait_all_children;
# Merge sub-files
my @all_cnt; 
for my $sfn ( @sub_fn ) {
        open F,'<',"$sfn.o" or &stopErr("[Err] Failed to open subfile [$sfn.o]\n");
	<F>; 
	my $i0 = -1; 
        while (<F>) {
		$i0 ++; 
		my $i = $i0 + $opts{'startColN'}; 
		chomp; 
		my @ta = split(/\t/, $_); 
		$ta[0] eq $ha[$i] or &stopErr("[Err] Different Individual ID [$i-$ta[0]] for $ha[$i] in file [$sfn.o]\nhs=[@ha]\n"); 
		$all_cnt[$i][0] //= $ta[0]; 
		$all_cnt[$i][0] eq $ta[0] or &stopErr("[Err] Different ID in all_cnt [$all_cnt[$i][0]] and [$ta[0]]\n"); 
		for ( my $j=1; $j<@ta; $j++ ) {
			$all_cnt[$i][$j] //= -1; 
			if ( $ta[$j] == -1 ) {
			} else {
				if ( $all_cnt[$i][$j] == -1 ) {
					$all_cnt[$i][$j] = $ta[$j]; 
				} else {
					$all_cnt[$i][$j] += $ta[$j]; 
				}
			}
		}
        }
        close F;
}
print STDOUT join("\t", qw/IndvID N_Num Typed_Num Het_Num Hom_Num Het_Ratio Hom_Ratio N_Ratio/)."\n"; 
for (my $i=$opts{'startColN'}; $i<@ha; $i++) {
	$ha[$i] eq $all_cnt[$i][0] or &stopErr("[Err] $ha[$i] VS. $all_cnt[$i][0] ?\n"); 
	my $tot        = $all_cnt[$i][1] + $all_cnt[$i][2]; 
	my $tot_typed  = $all_cnt[$i][2]; 
	my $rat_hete   = ($tot_typed > 0) ? ( $all_cnt[$i][3]/$tot_typed*100 ) : -1; 
	my $rat_homo   = ($tot_typed > 0) ? ( $all_cnt[$i][4]/$tot_typed*100 ) : -1; 
	my $rat_miss   = ($tot_typed > 0) ? ( $all_cnt[$i][1]/$tot      *100 ) : -1; 
	print STDOUT join("\t", @{$all_cnt[$i]}[0..4], $rat_hete, $rat_homo, $rat_miss)."\n"; 
}

# Delete temp_dir
&fileSunhh::_rmtree($wrk_dir);



sub type_tabGeno {
	if ( $_[0] eq './.' or $_[0] eq 'N/N' ) {
		return('miss'); 
	} elsif ( $_[0] =~ m!^([ATGC\*N]+)/([ATGC\*N]+)$! ) {
		if ( $1 eq $2 ) {
			return('homo'); 
		} else {
			return('hete'); 
		}
	} else {
		&stopErr("[Err] Unknown format of vcfTab genotype [$_[0]]\n"); 
	}
}# type_tabGeno
