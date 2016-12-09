#!/usr/bin/perl -w
# Cannot handle insertions. 
# For multi-cpu
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"noHeader!", 
	"cpuN:i", 
	"start_colN:i", 
	"skipInDel!", 
	"singleChar!", 
); 

$opts{'cpuN'} //= 20; 
$opts{'start_colN'} //= 3; 

my $help_txt = <<HH; 

perl $0 in.vcf.tab > in.vcf.cols

 -start_colN    [$opts{'start_colN'}]
 -cpuN          [$opts{'cpuN'}]
 -noHeader      [Boolean]

 -skipInDel     [Boolean] Skip lines with InDel alleles. 
 -singleChar    [Boolean] If given, 'AC'/'CA' in .cols format will be 'M'. 

# Mask indel genotypes to N, except single '*'. 
# Only accept bi-allele genotype for cols-format (single-character) genotype. 

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

my %dblist; 
my %b2dlist; 
{
	my @aa = (
	[qw/W A T/], 
	[qw/S C G/],
	[qw/M A C/],
	[qw/K G T/],
	[qw/R A G/],
	[qw/Y C T/]
	); 
	for my $tr (@aa) {
		my @bb = @$tr; 
		$dblist{$bb[0]} = [$bb[1], $bb[2]]; 
	}
	for my $tr (@aa) {
		$b2dlist{"$tr->[1]$tr->[2]"} = $tr->[0]; 
		$b2dlist{"$tr->[2]$tr->[1]"} = $tr->[0]; 
	}
	for my $tb (qw/A T G C/) {
		$b2dlist{$tb} = $tb; 
		$b2dlist{"$tb$tb"} = $tb; 
	}
}

my $fh = \*STDOUT; 
@ARGV and $fh = &openFH($ARGV[0]); 


unless ( defined $opts{'noHeader'} ) {
	my $l = <$fh>; 
	print STDOUT $l; 
}

my $pm = &LogInforSunhh::get_pm( $opts{'cpuN'} ); 

my %used; 
my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" ); 
for my $sfn ( @sub_fn ) {
	my $pid = $pm->start and next;
	open F,'<',"$sfn" or &stopErr("[Err] Failed to open subfile [$sfn]\n");
	open O,'>',"$sfn.o" or &stopErr("[Err] Failed to write subfile [$sfn.o]\n");
	while (<F>) {
		chomp;
		my @ta = &splitL("\t", $_); 
		my $is_bad_line = 0; 
		for my $tb ( @ta[$opts{'start_colN'} .. $#ta] ) {
			$tb = uc($tb); 
			if ($opts{'skipInDel'}) {
				$tb =~ m/[\*+]/ and do { $is_bad_line = 1; last; }; 
			}
			if ( $tb =~ s!^([ATGC\*])/\1$!$1! ) {
			} elsif ( $tb =~ s!^([ATGC\*])/([ATGC\*])$!$1$2! ) {
			} elsif ( $tb eq './.' or $tb eq 'N/N' ) {
				$tb = 'N'; 
			} elsif ( $tb =~ m!^[ATGC\*]+/[ATGC\*]+$! ) {
				if ( $opts{'skipInDel'} ) {
					$is_bad_line = 1; 
					last; 
				}
				unless ( defined $used{'skip_geno'}{$tb} ) {
					&tsmsg("[Wrn] Mask genotype [$tb]\n"); 
					$used{'skip_geno'}{$tb} = 1; 
				}
				$tb = 'N'; 
			} elsif ( $tb =~ m!^[ATGCN\*]$! ) {
				; 
			} elsif ( defined $dblist{$tb} ) {
				$tb = $dblist{$tb}[0] . $dblist{$tb}[1]; 
			} else {
				unless ( defined $used{'skip_geno'}{$tb} ) {
					&tsmsg("[Wrn] Mask genotype [$tb]\n"); 
					$used{'skip_geno'}{$tb} = 1; 
				}
				$tb = 'N'; 
			}
			$tb eq '*' and $tb = '-'; 
			if ($opts{'singleChar'}) {
				unless ( $tb eq '-' or $tb eq 'N' ) {
					defined $b2dlist{$tb} or &stopErr("[Err] Unknown genotype [$tb]\n"); 
					$tb = $b2dlist{$tb}; 
				}
			}
		}
		$is_bad_line == 0 or next; 
		print O join("\t", @ta)."\n"; 
	}
	close O;
	close F;
	$pm->finish;
}
$pm->wait_all_children;

for my $sfn ( @sub_fn ) {
	open F,'<',"$sfn.o" or &stopErr("[Err] Failed to open subfile [$sfn.o]\n");
	while (<F>) {
		print STDOUT $_;
	}
	close F;
}
&fileSunhh::_rmtree($wrk_dir); 


