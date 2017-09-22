#!/usr/bin/perl 
# 2016-11-27 This is too too slow for big data, so I want to add multi-threading for it. 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"startColN:i", # 2
	"indel_asN!", # Rule 1 : Treat all indel as N missing; 
	"hete_asN!",  # Rule 2 : Treat heterozygous as missing; 
	"maxMiss:f",  # >= 0.05 
	"only_count!", # Do not do the filtering if given. 
	"cpuN:i", # 
); 

$opts{'maxMiss'} //= 0.05; 
$opts{'startColN'} //= 2; 
$opts{'cpuN'}  //= 1; 
$opts{'cpuN'} < 1 and $opts{'cpuN'} = 1; 


my $maxAllowMiss = $opts{'maxMiss'} * 100; 

sub usage {
	print STDERR <<HH; 

perl $0 in.snp > in.snp.filtered

-help 
-startColN        [$opts{'startColN'}]
-maxMiss          [$opts{'maxMiss'}]  missing rate should be <= maxMiss 
-indel_asN        [Bool] If given, genotype with '*|+' will be replaced with N
-hete_asN         [Bool] If given, genotype with m/^[ATGCN]{2,}$/ will be replaced with N 
-only_count       [Bool] Only count the missing rate instead of filter the table. 

-cpuN             [1] If >= 2, I'll split the input file into \$cpuN files to process. 

HH
	exit(1); 
}

-t and !@ARGV and &usage(); 
$opts{'help'} and &usage(); 

if ( $opts{'cpuN'} > 1 ) {
	use Parallel::ForkManager; 
	my $MAX_PROCESSES = $opts{'cpuN'}; 
	my $pm = new Parallel::ForkManager($MAX_PROCESSES);
	my $all_lineN = 0; 
	my $wrk_dir = &fileSunhh::new_tmp_dir(); 
	mkdir($wrk_dir) or &stopErr("[Err] Failed to creat dir [$wrk_dir]\n"); 
	my $ofh = &openFH("$wrk_dir/single_input", '>'); 
	&tsmsg("[Msg] Creating single input file [$wrk_dir/single_input]\n"); 
	while (<>) {
		$all_lineN ++; 
		print {$ofh} $_; 
	}
	close($ofh); 
	my $per_lineN = $all_lineN / $opts{'cpuN'}; 
	int($per_lineN) == $per_lineN or $per_lineN = int($per_lineN) + 1; 
	open F,'<', "$wrk_dir/single_input" or &stopErr("[Err] Failed to open file [$wrk_dir/single_input]\n"); 
	my $end_lineN = $per_lineN; 
	my $o_id = int( $end_lineN / $per_lineN ); 
	open O,'>', "$wrk_dir/out_${o_id}" or &stopErr("[Err] Failed to create dvd file [$wrk_dir/out_${o_id}]\n"); 
	while (<F>) {
		unless ( $. <= $end_lineN ) {
			$end_lineN += $per_lineN; 
			$. <= $end_lineN or &stopErr("[Err] \$. = $. and end_lineN = $end_lineN\n"); 
			$o_id = int( $end_lineN / $per_lineN ); 
			close(O); 
			open O,'>', "$wrk_dir/out_${o_id}" or &stopErr("[Err] Failed to create dvd file [$wrk_dir/out_${o_id}]\n"); 
			&tsmsg("[Msg] Writing dvd file [$wrk_dir/out_${o_id}]\n"); 
		}
		print O $_; 
	}
	close F; 
	close O; 
	for (my $i=1; $i<=$o_id; $i++) {
		my $pid = $pm->start and next; 
		&tsmsg("[Msg]   Processing dvd file [$wrk_dir/out_${i}]\n"); 
		open I1, '<', "$wrk_dir/out_${i}" or &stopErr("[Err] Failed to read in dvd file [$wrk_dir/out_${i}]\n"); 
		open O1, '>', "$wrk_dir/out_${i}_final" or &stopErr("[Err] Failed to creat dvd out [$wrk_dir/out_${i}_final]\n"); 
		while (<I1>) {
			&process_in_line($_, \*O1); 
		}
		close O1; 
		close I1; 
		&tsmsg("[Msg]   Finished dvd file [$wrk_dir/out_${i}]\n"); 
		$pm->finish; 
	}
	$pm->wait_all_children; 
	for (my $i=1; $i<=$o_id; $i++) {
		open I2,'<',"$wrk_dir/out_${i}_final" or die; 
		while (<I2>) {
			print STDOUT $_; 
		}
		close I2; 
	}
	&fileSunhh::_rmtree($wrk_dir); 
} else {
	while (<>) {
		&process_in_line($_, \*STDOUT); 
	}
}

sub process_in_line {
	$_ = $_[0]; 
	s/[^\t\S]+$//; 
	
	my @ta = split(/\t/, $_); 
	my ($chr, $pos) = @ta[0,1]; 
	if ($chr eq 'chr') {
		if ( $opts{'only_count'} ) {
			print {$_[1]} join("\t", qw/chr pos NmissRate/)."\n"; 
		} else {
			print {$_[1]} "$_\n"; 
		}
		return; 
	}
	
	# Counting 
	my $missingCnt = 0; 
	my $totalCnt = 0; 
	for my $tb (@ta[$opts{'startColN'}..$#ta]) {
		$tb = uc($tb); 
		$opts{'indel_asN'} and $tb =~ m/[*+]/ and $tb = 'N';  # R1 
		$opts{'hete_asN'} and $tb =~ m/^[ATGCN*]{2,}$/ and $tb = 'N'; # R2 
		$tb eq 'N' and $missingCnt++; 
		$totalCnt ++; 
	}
	my $missingRate = int($missingCnt/$totalCnt*10000+0.5)/100; 
	
	# Filtering
	if ( $opts{'only_count'} ) {
		print {$_[1]} join("\t", $chr, $pos, $missingRate)."\n"; 
	} else {
		if ( $missingRate <= $maxAllowMiss) {
			print {$_[1]} "$_\n"; 
		}
	}
}

