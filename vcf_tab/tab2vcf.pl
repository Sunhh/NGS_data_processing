#!/usr/bin/perl
use strict; 
use warnings; 
use Getopt::Long; 
use LogInforSunhh; 
use fileSunhh; 
use fastaSunhh; 
use SNP_tbl; 
my $fas_obj = fastaSunhh->new(); 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"ref_fas:s", # db/ITAG2.3_genomic.fa 
	"in_tab:s", # 
	"cpuN:i", 
); 

$opts{'cpuN'} //= 1; 
my $pm; 
$opts{'cpuN'} > 1 and $pm = &LogInforSunhh::get_pm( $opts{'cpuN'} );


my $help_txt = <<HH; 

perl $0 -ref_fas db/ITAG2.3_genomic.fa   -in_tab Slycopersicoides_snps_rmfilt02V.tab

-help      
-cpuN    [$opts{'cpuN'}]

HH

( defined $opts{'ref_fas'} and defined $opts{'in_tab'} ) or &LogInforSunhh::usage($help_txt); 

my %seq = %{ $fas_obj->save_seq_to_hash( 'faFile'=>$opts{'ref_fas'} ) }; 
for (keys %seq) { $seq{$_}{'seq'} =~ s!\s!!g; $seq{$_}{'len'} = length($seq{$_}{'seq'}); }

my $o_header = <<OH; 
##fileformat=VCFv4.1
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="The phred-scaled genotype likelihoods rounded to the closest integer">
OH
for my $k1 (sort keys %seq) {
	my $ltxt = "##contig=<ID=$k1,length=$seq{$k1}{'len'}>\n"; 
	$o_header .= $ltxt; 
}
{
	my $absP = &fileSunhh::_abs_path($opts{'in_tab'}); 
	$o_header .= "##reference=file://$absP\n"; 
}


print STDOUT "$o_header"; 
my $fh = &openFH( $opts{'in_tab'}, '<' ); 
my $has_header = 0; 
{
	$_ = <$fh>; 
	chomp; 
	my @ta = split(/\t/, $_); 
	if ( $ta[1] =~ m/^pos$/i ) {
		if ( $has_header == 0 ) {
			$has_header = 1; 
			print STDOUT join("\t", "#CHROM", "POS", qw/ID REF ALT QUAL FILTER INFO FORMAT/, @ta[3..$#ta])."\n"; 
		} else {
			&stopErr("[Err] repeat header.\n$_\n"); 
		}
	}
}
if ( $opts{'cpuN'} <= 1 ) {
	while (my $il = <$fh>) {
		my ($ol) = &get_line($il); 
		defined $ol or next; 
		print STDOUT $ol; 
	}
	close($fh); 
} else {
	my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 );
	my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" ); 
	close($fh); 
	for my $sfn (@sub_fn) {
		my $pid = $pm->start and next;
		open F,'<',"$sfn" or die;
		open O,'>',"$sfn.o" or die;
		while (<F>) {
			my ($ol) = &get_line($_); 
			defined $ol or next; 
			print O $ol; 
		}
		close O; 
		close F; 
		$pm->finish;
	}
	$pm->wait_all_children;
	for my $sfn (@sub_fn) {
		open F,'<',"$sfn.o" or die; 
		while ( <F> ) { print STDOUT $_; }
		close F; 
	}
	&fileSunhh::_rmtree($wrk_dir); 
}

sub get_line {
	my $il = $_[0]; 
	chomp($il); 
	my @ta = split(/\t/, $il); 
	my $mrkID = "$ta[0]_$ta[1]"; 
	my $refBase = uc( substr( $seq{$ta[0]}{'seq'}, $ta[1]-1, 1 ) ); 
	my %alleles; 
	for my $tb (@ta[3 .. $#ta]) {
		my @idv_al = &SNP_tbl::tab_allele($tb); 
		$idv_al[0][0] eq '.' and next; 
		for my $tc (@idv_al) {
			$alleles{ $tc->[0] } ++; 
		}
	}
	delete( $alleles{$refBase} ); 
	my @sort_al = sort { $alleles{$b} <=> $alleles{$a} || $a cmp $b } keys %alleles; 
	my @arr_ALT; 
	my %geno2num; 
	$geno2num{$refBase} = 0; 
	for (my $i=0; $i<@sort_al; $i++) {
		$geno2num{$sort_al[$i]} = $i+1; 
		push( @arr_ALT, $sort_al[$i] ); 
	}
	my $txt_ALT = ( scalar(@arr_ALT) > 0 ) ? join(',', @arr_ALT) : '.' ; 
	for my $tb ( @ta[3 .. $#ta] ) {
		my @idv_al = &SNP_tbl::tab_allele($tb); 
		$idv_al[0][0] eq '.' and do { $tb = './.'; next; }; 
		$#idv_al == 0 and push(@idv_al, $idv_al[0]); 
		$tb = join('/', sort { $a <=> $b ; } map { $geno2num{$_->[0]} } @idv_al ); 
	}
	return( join("\t", $ta[0], $ta[1], '.', $refBase, $txt_ALT, '.', '.', '.', 'GT', @ta[3 .. $#ta])."\n" ); 
}# get_line() 

