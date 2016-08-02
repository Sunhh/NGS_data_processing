#!/usr/bin/perl
use strict; 
use warnings; 
use Getopt::Long; 
use fileSunhh; 
use fastaSunhh; 
use SNP_tbl; 
my $fas_obj = fastaSunhh->new(); 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"ref_fas:s", # db/ITAG2.3_genomic.fa 
	"in_tab:s", # 
); 


my $help_txt = <<HH; 

perl $0 -ref_fas db/ITAG2.3_genomic.fa   -in_tab Slycopersicoides_snps_rmfilt02V.tab

-help      

HH

( defined $opts{'ref_fas'} and defined $opts{'in_tab'} ) or &LogInforSunhh::usage($help_txt); 

my %seq = %{ $fas_obj->save_seq_to_hash( 'faFile'=>$opts{'ref_fas'} ) }; 
for (keys %seq) { $seq{$_}{'seq'} =~ s!\s!!g; $seq{$_}{'len'} = length($seq{$_}{'seq'}); }

my $o_header = <<OH; 
##fileformat=VCFv4.1
###ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
OH

print "$o_header"; 
my $fh = &openFH( $opts{'in_tab'}, '<' ); 
my $has_header = 0; 
while (<$fh>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ( $ta[1] =~ m/^pos$/i ) {
		if ( $has_header == 0 ) {
			$has_header = 1; 
			print join("\t", "#CHROM", "POS", qw/ID REF ALT QUAL FILTER INFO FORMAT/, @ta[2..$#ta])."\n"; 
			next; 
		} else {
			&stopErr("[Err] repeat header.\n$_\n"); 
		}
	}
	my $mrkID = "$ta[0]_$ta[1]"; 
	my $refBase = uc( substr( $seq{$ta[0]}{'seq'}, $ta[1]-1, 1 ) ); 
	my %alleles; 
	for my $tb (@ta[2 .. $#ta]) {
		my @idv_al = &SNP_tbl::tab_allele($tb); 
		$idv_al[0][0] eq '.' and next; 
		for my $tc (@idv_al) {
			$alleles{ $tc->[0] } ++; 
		}
	}
	delete( $alleles{$refBase} ); 
	my @sort_al = sort { $alleles{$b} <=> $alleles{$a} } keys %alleles; 
	my @arr_ALT; 
	my %geno2num; 
	$geno2num{$refBase} = 0; 
	for (my $i=0; $i<@sort_al; $i++) {
		$geno2num{$sort_al[$i]} = $i+1; 
		push( @arr_ALT, $sort_al[$i] ); 
	}
	my $txt_ALT = ( scalar(@arr_ALT) > 0 ) ? join(',', @arr_ALT) : '.' ; 
	for my $tb ( @ta[2 .. $#ta] ) {
		my @idv_al = &SNP_tbl::tab_allele($tb); 
		$idv_al[0][0] eq '.' and do { $tb = './.'; next; }; 
		$#idv_al == 0 and push(@idv_al, $idv_al[0]); 
		$tb = join('/', sort { $a <=> $b ; } map { $geno2num{$_->[0]} } @idv_al ); 
	}
	print join("\t", $ta[0], $ta[1], '.', $refBase, $txt_ALT, '.', '.', '.', 'GT', @ta[2 .. $#ta])."\n"; 
}
close($fh); 

