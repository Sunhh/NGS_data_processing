#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"startColN:i", # 2 
	"showTime:i", 
	"cpuN:i", # 20
); 
use SNP_tbl; 
my $st_obj = SNP_tbl->new(); 

$opts{'startColN'} //= 2; 
$opts{'cpuN'}      //= 1; 

my $help_txt = <<HH ; 

perl $0 in_snp.tbl > in_snp.fasta

Need a headere line. 

-help
-startColN       [$opts{'startColN'}]

-showTime        [0] Will report time when 'showTime' number of lines have been processed. 

-cpuN            [$opts{'cpuN'}]

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 
$opts{'showTime'} //= 0; 

my $pm; 
$pm = &LogInforSunhh::get_pm( $opts{'cpuN'} ); 

my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 );
my $fh = \*STDIN; 
@ARGV and $fh = &openFH( $ARGV[0] , '<' ); 

my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 1, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" );
close($fh); 

for my $sfn (@sub_fn) {
	my $pid = $pm->start and next;
	&tbl2seq( $sfn, "$sfn.o" ); 
	$pm->finish;
}
$pm->wait_all_children;

my (@okey, %oseq); 
for my $sfn (@sub_fn) {
	my %sss = %{ $fs_obj->save_seq_to_hash('faFile'=>"$sfn.o") }; 
	@okey == 0 and @okey = sort { $sss{$a}{'Order'} <=> $sss{$b}{'Order'} } keys %sss; 
	for my $tk (@okey) {
		$sss{$tk}{'seq'} =~ s/\s//g; 
		$oseq{$tk} .= $sss{$tk}{'seq'}; 
	}
}
for my $tk (@okey) {
	$oseq{$tk} =~ s!\s!!g; 
	$oseq{$tk} =~ s!(\S{100})!$1\n!g; 
	chomp($oseq{$tk}); 
	print STDOUT ">$tk\n$oseq{$tk}\n"; 
}
&fileSunhh::_rmtree($wrk_dir); 

sub tbl2seq {
	my $fn = shift; 
	my $ofn = shift; 
	my $ifh = &openFH($fn, '<'); 
	my $hl = <$ifh>; 
	chomp($hl); 
	my @hh = split(/\t/, $hl); 
	my @seq; 
	while (<$ifh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		for (my $i=$opts{'startColN'}; $i<@ta; $i++) {
			if ( $ta[$i] =~ m!/! ) {
				if ( $ta[$i] =~ s!^(\S+)/\1$!$1! ) {
				} elsif ( $ta[$i] =~ s!^([ATGC])/([ATGC])$!$1$2! ) {
				} else {
					$ta[$i] = "N"; 
				}
			}
			$ta[$i] eq '.' and $ta[$i] = "N"; 
			if ( $ta[$i] eq '*' ) {
				$ta[$i] = '-'; 
			} else {
				$ta[$i] = $st_obj->SingleChar($ta[$i], 'maxAlleleN'=>2); 
			}
			$seq[$i] .= $ta[$i]; 
		}
	}
	close ($ifh); 

	my $ofh = &openFH($ofn, '>'); 
	for (my $i=$opts{'startColN'}; $i<@hh; $i++) {
		print {$ofh} ">$hh[$i]\n"; 
		$seq[$i] =~ s/\s//g; 
		$seq[$i] =~ s/(\S{100})/$1\n/g; 
		chomp($seq[$i]); 
		print ${ofh} "$seq[$i]\n"; 
	}
	close($ofh); 
}





