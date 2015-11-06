#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use SNP_tbl; 
my $st_obj = SNP_tbl->new(); 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"minMAF:f", 
	"maxMissR:f", 
	"only2allele!", 
	"startColN:i", 
	"col_list:s", 
); 

$opts{'startColN'} //= 2; 
my $geno_col = $opts{'startColN'}; 
$opts{'minMAF'} //= 0.05; # OK >= minMAF; 
$opts{'maxMissR'} //= 0.40; # OK <= maxMissR; 

!@ARGV and !(defined $opts{'col_list'}) and &usage(); 
$opts{'help'} and &usage(); 

sub usage {
print STDERR <<HH; 

perl $0 in.cols 
# Will write a new file called 'in.cols.filter' 

-col_list    [''] list of .cols files. Will output files with '.filter' suffix; 

-help
-minMAF      [$opts{'minMAF'}]
-maxMissR    [$opts{'maxMissR'}]
-only2allele Should be given. 
-startColN   [$opts{'startColN'}]

HH
	exit(1); 
}

if (defined $opts{'col_list'}) {
	@ARGV = (); 
	my $l_fh = &openFH($opts{'col_list'}, '<'); 
	while (<$l_fh>) {
		m/^\s*$/ and next; 
		m/^\s*#/ and next; 
		chomp; 
		my @ta = split(/\t/, $_); 
		push(@ARGV, $ta[0]); 
	}
	close ($l_fh); 
}

for my $inF (@ARGV) {

&tsmsg("[Rec] Processing [$inF]\n"); 
my $in_fh = &openFH($inF, '<'); 
my $out_fh = &openFH("$inF.filter", '>'); 

my $header = <$in_fh>; 
print {$out_fh} "$header"; 
my (%al_cnt, $al_tot, $miss_cnt); 
my @al_base; 
LINE: 
while (<$in_fh>) {
	$. % 1e6 == 1 and &tsmsg("[Msg]   Loading $. line.\n"); 
	chomp; 
	my @ta = split(/\t/, $_); 
	%al_cnt = (); $al_tot = 0; $miss_cnt = 0; 
	for my $tb (@ta[$geno_col .. $#ta]) {
		$tb eq '-' and $tb = 'N'; 
		$tb eq 'N' and do { $miss_cnt ++; next; }; 
		# $tb = ($opts{'SNP_tbl::dna_d2b($tb)'}) ? $st_obj->SingleChar( $tb, 'maxAlleleN'=>2, 'onlyATGC'=>0 ) : $st_obj->SingleChar( $tb, 'maxAlleleN'=>0, 'onlyATGC'=>0 ); 
		$tb = $st_obj->SingleChar( $tb, 'maxAlleleN'=>2, 'onlyATGC'=>0 ) ; 
		$tb eq 'N' and next LINE; # Contain abnormal genotype. 
		@al_base = &SNP_tbl::dna_d2b($tb); 
		$al_base[1] //= $al_base[0]; 
		$al_cnt{ $al_base[0] } ++; 
		$al_cnt{ $al_base[1] } ++; 
		$al_tot += 2; 
	}
	my @tk = sort { $al_cnt{$b} <=> $al_cnt{$a} || $a cmp $b } keys %al_cnt; 
	# defined $tk[1] or next; 
	if ( $opts{'only2allele'} ) {
		scalar(@tk) == 2 or next; 
	}
	$miss_cnt <= ($#ta-$geno_col+1) * $opts{'maxMissR'} or next LINE; 
	my $maf = ( defined $tk[1] ) ? $al_cnt{$tk[1]}/$al_tot : 0 ; 
	$maf >= $opts{'minMAF'} or next LINE; 
	
	print {$out_fh} join("\t", @ta)."\n"; 
}
close($out_fh); 
close($in_fh); 
}

&tsmsg("[Rec] All done.\n"); 
