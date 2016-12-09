#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use SNP_tbl; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"startColN:i", # 2 
	"noHeader!", 
	"cpuN:i", # 20 
); 
$opts{'startColN'} //= 2; 
$opts{'cpuN'} //= 20; 
$opts{'cpuN'} < 1 and $opts{'cpuN'} = 1; 

my $pm = &LogInforSunhh::get_pm( $opts{'cpuN'} ); 
my %d2b_list = &SNP_tbl::get_diploid_d2b(); 

my $help_txt = <<HH; 

Count N (missing), heterozygous, and homozygous sites per individual. 
The output format is : qw/IndvID N_Num Typed_Num Het_Num Hom_Num Het_Ratio Hom_Ratio N_Ratio/

perl $0 in.snp > in.snp.cntNHH

 -help               

 -cpuN               [$opts{'cpuN'}]
 -startColN          [$opts{'startColN'}]
 -noHeader           [Bool]

Genotype column start from colN=$opts{'startColN'}
Do not parse the first line. 

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my $fh = \*STDOUT; 
@ARGV > 0 and $fh = &openFH( $ARGV[0], '<' ); 

my @ha; 
unless ($opts{'noHeader'}) {
	my $head = <>; 
	chomp($head); 
	@ha=split(/\t/, $head); 
}

my @cnt_N; 
my @cnt_hete; 
my @cnt_homo; 

while (<>) { 
	($.-1) % 1e6 == 1 and &tsmsg("[Msg] $. lines.\n"); 
	chomp; 
	my @ta = split(/\t/, $_); 
} 
print STDOUT join("\t", qw/IndvID N_Num Typed_Num Het_Num Hom_Num Het_Ratio Hom_Ratio N_Ratio/)."\n"; 
for (my $i=$opts{'startColN'}; $i<@ha; $i++) { 
	$cnt_N[$i] //= 0; 
	$cnt_hete[$i] //= 0; 
	$cnt_homo[$i] //= 0; 
	my $tot = $cnt_N[$i]+$cnt_hete[$i]+$cnt_homo[$i]; 
	my $tot_typed = $cnt_hete[$i]+$cnt_homo[$i]; 
	my $rat_hete = ($tot_typed > 0) ? sprintf("%02.02f", $cnt_hete[$i]/$tot_typed*100) : -1 ; 
	my $rat_homo = ($tot_typed > 0) ? sprintf("%02.02f", $cnt_homo[$i]/$tot_typed*100) : -1 ; 
	my $rat_miss = ($tot       > 0) ? sprintf("%02.02f", $cnt_N[$i]   /$tot      *100) : -1 ; 
	print "$ha[$i]\t$cnt_N[$i]\t$tot_typed\t$cnt_hete[$i]\t$cnt_homo[$i]\t$rat_hete\t$rat_homo\t$rat_miss\n";
}


# separate files
my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 );
my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" );
# Process sub-files
for my $sfn ( @sub_fn ) {
        my $pid = $pm->start and next;
        open F,'<',"$sfn" or &stopErr("[Err] Failed to open subfile [$sfn]\n");
        open O,'>',"$sfn.o" or &stopErr("[Err] Failed to write subfile [$sfn.o]\n");
        while (<F>) {
                chomp;
                my @ta = &splitL("\t", $_);
		if ($opts{'noHeader'} and @ha == 0) {
			for (my $i=0; $i<@ta; $i++) {
				$ha[$i] = "Col_$i"; 
			}
		}
		for (my $i=$opts{'startColN'}; $i<@ta; $i++) { 
			($ta[$i] eq 'N' or $ta[$i] eq 'n') and do { $cnt_N[$i]++; next; }; 
			( $ta[$i] =~ m/^[ATGC*]$/i or $ta[$i] =~ m/\+/ ) and do { $cnt_homo[$i]++; next; }; # The '*' and sites with \+ are treated as homozygous. I don't like genotype like 'A*', so I want to remove it before calculation. 
			(&SNP_tbl::dna_d2b( &SNP_tbl::dna_b2d($ta[$i]) )) > 1 and do { $cnt_hete[$i]++; next; }; # Heterozygous. 
			&tsmsg("[Wrn] Weired genotype [$ta[$i]] is treated as homozygous.\n"); 
			$cnt_homo[$i] ++; 
			# $ta[$i]=~m/^[ATGC*]$|^[ATGC]\+[ATGC]+$/ and $cnt[$i]++; 
		} 
                print O join("\t", map { ( $_ eq 'NA' ) ? '' : $ta[$_] ; } @cur_ColNs)."\n";
        }
        close O;
        close F;
        $pm->finish;
}
sub cnt_stat {
}
$pm->wait_all_children;
# Merge sub-files
print STDOUT $header_txt;
for my $sfn ( @sub_fn ) {
        open F,'<',"$sfn.o" or &stopErr("[Err] Failed to open subfile [$sfn.o]\n");
        while (<F>) {
                print STDOUT $_;
        }
        close F;
}
# Delete temp_dir
&fileSunhh::_rmtree($wrk_dir);



