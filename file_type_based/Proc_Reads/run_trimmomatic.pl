#!/usr/bin/perl
use strict; 
use FindBin; (my $REPO = $FindBin::RealBin) =~ s{(/NGS_data_processing)(/.*)?$}{$1};  # portable repo root
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"fq1_pref:s@", "fq1_suff:s", "fq1_pref_list:s", 
	"fq2_pref:s@", "fq2_suff:s", "fq2_pref_list:s", 
	"jar_trimmo:s", 
	"adp_PE:s", "adp_SE:s", "adp_polyAT:s", 
	"paraSE_trimmo:s", "paraPE_trimmo:s", "paraJar_trimmo:s", "paraSE_trimmo_polyAT:s", 
	"exe_java:s", 
	"trimPolyAT!", 
	"bgzip!", "exe_bgzip:s", 
); 

$opts{'fq1_suff'} //= '_R1.ndupB'; 
$opts{'fq2_suff'} //= '_R2.ndupB'; 
$opts{'jar_trimmo'} //= "$REPO/file_type_based/Proc_Reads/trimmomatic/trimmomatic.jar"; 
$opts{'paraJar_trimmo'} //= '-threads 20'; 
$opts{'exe_java'} //= 'java'; 
$opts{'exe_bgzip'} //= 'bgzip'; 

# Adapter fasta files: default to the in-package copies; override with -adp_PE / -adp_SE / -adp_polyAT .
$opts{'adp_PE'}     //= "$REPO/file_type_based/Proc_Reads/trimmomatic/adapters/TruSeq3-PE-2.fa"; 
$opts{'adp_SE'}     //= "$REPO/file_type_based/Proc_Reads/trimmomatic/adapters/TruSeq3-SE.fa"; 
$opts{'adp_polyAT'} //= "$REPO/file_type_based/Proc_Reads/polyAT_adp.fa"; 

# The default parameter combinations are the historical ones (kept identical to the retired
# run_trimmoPE.pl / run_trimmoSE.pl), with only the adapter path taken from the -adp_* options.
# Override the whole string with -paraPE_trimmo / -paraSE_trimmo if needed.
$opts{'paraPE_trimmo'}        //= "ILLUMINACLIP:$opts{'adp_PE'}:2:30:10:1:TRUE SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40"; 
$opts{'paraSE_trimmo'}        //= "ILLUMINACLIP:$opts{'adp_SE'}:2:30:10 SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40"; 
$opts{'paraSE_trimmo_polyAT'} //= "ILLUMINACLIP:$opts{'adp_polyAT'}:2:30:7 MINLEN:40"; 

if ( $opts{'trimPolyAT'} ) {
	$opts{'paraSE_trimmo'} = $opts{'paraSE_trimmo_polyAT'}; 
}

defined $opts{'fq1_pref_list'} and push( @{$opts{'fq1_pref'}}, @{ &load_prefList( $opts{'fq1_pref_list'} ) } ); 
defined $opts{'fq2_pref_list'} and push( @{$opts{'fq2_pref'}}, @{ &load_prefList( $opts{'fq2_pref_list'} ) } ); 
$opts{'fq1_pref'} //= []; 
$opts{'fq2_pref'} //= []; 
my @fq1_pref = @{ $opts{'fq1_pref'} }; 
my @fq2_pref = @{ $opts{'fq2_pref'} }; 

my $help_txt = <<HH; 
perl $0 -fq1_pref inFq1_pref -fq1_suff _R1.ndupB   [ -fq2_pref inFq2_pref -fq2_suff _R2.ndupB ]

For each index i, if both fq1_pref[i] and fq2_pref[i] are given -> Trimmomatic PE, else SE.
Input file = <pref><suff>. PE outputs: <fq1_pref>_pTr_R1.fq <fq1_pref>_sTr_R1.fq <fq2_pref>_pTr_R2.fq <fq2_pref>_sTr_R2.fq ;
SE output: <pref>_sTr_R1.fq (fq1) or <pref>_sTr_R2.fq (fq2). (Use fq1_pref == fq2_pref to share one prefix.)

-fq1_pref_list          [filename] one prefix per line (col1); appended to -fq1_pref
-fq2_pref_list          [filename] same, for read-2 prefixes

-jar_trimmo             [$opts{'jar_trimmo'}]
-exe_java               [$opts{'exe_java'}]
-paraJar_trimmo         [$opts{'paraJar_trimmo'}]

Adapter fasta (in-package by default; parameter-specifiable):
-adp_PE                 [$opts{'adp_PE'}]
-adp_SE                 [$opts{'adp_SE'}]
-adp_polyAT             [$opts{'adp_polyAT'}]

Full Trimmomatic step strings (default = historical combos, adapter from -adp_*):
-paraPE_trimmo          [$opts{'paraPE_trimmo'}]
-paraSE_trimmo          [$opts{'paraSE_trimmo'}]
-trimPolyAT             [Boolean] use the poly-A/T SE string instead of -paraSE_trimmo
-paraSE_trimmo_polyAT   [$opts{'paraSE_trimmo_polyAT'}]

Compression:
-bgzip                  [Boolean] bgzip the output .fq files after each run
-exe_bgzip              [$opts{'exe_bgzip'}]

HH
$opts{'help'} and &LogInforSunhh::usage($help_txt); 
( @fq1_pref > 0 or @fq2_pref > 0 ) or &LogInforSunhh::usage($help_txt); 


for (my $i=0; $i<@fq1_pref; $i++) {
	$fq1_pref[$i] //= ''; 
	$fq2_pref[$i] //= ''; 
	my @outFq; 
	if ( $fq1_pref[$i] ne '' and $fq2_pref[$i] ne '') {
		my $fq1 = "$fq1_pref[$i]" . "$opts{'fq1_suff'}"; 
		my $fq2 = "$fq2_pref[$i]" . "$opts{'fq2_suff'}"; 
		@outFq = ( "$fq1_pref[$i]_pTr_R1.fq", "$fq1_pref[$i]_sTr_R1.fq", "$fq2_pref[$i]_pTr_R2.fq", "$fq2_pref[$i]_sTr_R2.fq" ); 
		&exeCmd_1cmd( "$opts{'exe_java'} -jar $opts{'jar_trimmo'} PE $opts{'paraJar_trimmo'}   $fq1   $fq2   $outFq[0] $outFq[1] $outFq[2] $outFq[3] $opts{'paraPE_trimmo'}" ); 
	} elsif ( $fq1_pref[$i] ne '' ) {
		@outFq = ( "$fq1_pref[$i]_sTr_R1.fq" ); 
		&exeCmd_1cmd( "$opts{'exe_java'} -jar $opts{'jar_trimmo'} SE $opts{'paraJar_trimmo'}   $fq1_pref[$i]$opts{'fq1_suff'}   $outFq[0]   $opts{'paraSE_trimmo'}" ); 
	} elsif ( $fq2_pref[$i] ne '' ) {
		@outFq = ( "$fq2_pref[$i]_sTr_R2.fq" ); 
		&exeCmd_1cmd( "$opts{'exe_java'} -jar $opts{'jar_trimmo'} SE $opts{'paraJar_trimmo'}   $fq2_pref[$i]$opts{'fq2_suff'}   $outFq[0]   $opts{'paraSE_trimmo'}" ); 
	} else {
		&tsmsg("[Wrn] Skip null line [$i].\n"); 
		next; 
	}
	if ( $opts{'bgzip'} ) {
		for my $ofq ( @outFq ) {
			-e $ofq or do { &tsmsg("[Wrn] Output [$ofq] not found; skip bgzip.\n"); next; }; 
			&exeCmd_1cmd( "$opts{'exe_bgzip'} -f $ofq" ); 
		}
	}
}

sub load_prefList {
	my ($pref_list) = @_; 
	my @back; 
	my $fh = &openFH($pref_list, '<'); 
	while ( &wantLineC($fh) ) {
		my @ta = &splitL("\t", $_); 
		# $ta[0] eq '' and next; 
		push(@back, $ta[0]); 
	}
	close($fh); 
	return(\@back); 
}

