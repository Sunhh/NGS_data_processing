#!/usr/bin/perl
use strict; 
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
	"paraSE_trimmo:s", "paraPE_trimmo:s", "paraJar_trimmo:s", 
	"exe_java:s", 
); 

$opts{'fq1_suff'} //= '_R1.ndupB'; 
$opts{'fq2_suff'} //= '_R2.ndupB'; 
$opts{'jar_trimmo'} //= '/home/Sunhh/src/Assemble/Trimmomatic/Trimmomatic-0.32/trimmomatic-0.32.jar'; 
$opts{'paraJar_trimmo'} //= '-threads 20'; 
$opts{'paraSE_trimmo'} //= 'ILLUMINACLIP:/home/Sunhh/src/Assemble/Trimmomatic/Trimmomatic-0.32/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40';
$opts{'paraPE_trimmo'} //= 'ILLUMINACLIP:/home/Sunhh/src/Assemble/Trimmomatic/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40'; 
$opts{'exe_java'} //= 'java'; 

defined $opts{'fq1_pref_list'} and push( @{$opts{'fq1_pref'}}, @{ &load_prefList( $opts{'fq1_pref_list'} ) } ); 
defined $opts{'fq2_pref_list'} and push( @{$opts{'fq2_pref'}}, @{ &load_prefList( $opts{'fq2_pref_list'} ) } ); 
$opts{'fq1_pref_list'} //= []; 
$opts{'fq2_pref_list'} //= []; 
my @fq1_pref = @{ $opts{'fq1_pref_list'} }; 
my @fq2_pref = @{ $opts{'fq2_pref_list'} }; 

my $help_txt = <<HH; 
perl $0 -fq1_pref inFq1_pref -fq1_suff _R1.ndupB   [ -fq2_pref inFq2_pref -fq2_suff _R2.ndupB ]

-fq1_pref_list
-fq2_pref_list
-jar_trimmo             [$opts{'jar_trimmo'}]
-paraSE_trimmo          [$opts{'paraSE_trimmo'}]
-paraPE_trimmo          [$opts{'paraPE_trimmo'}]
-paraJar_trimmo         [$opts{'paraJar_trimmo'}]
-exe_java               [$opts{'exe_java'}]

HH
$opts{'help'} and &LogInforSunhh::usage($help_txt); 
( @fq1_pref > 0 or @fq2_pref > 0 ) or &LogInforSunhh::usage($help_txt); 


for (my $i=0; $i<@fq1_pref; $i++) {
	$fq1_pref[$i] //= ''; 
	$fq2_pref[$i] //= ''; 
	if ( $fq1_pref[$i] ne '' and $fq2_pref[$i] ne '') {
		my $fq1 = "$fq1_pref[$i]" . "$opts{'fq1_suff'}"; 
		my $fq2 = "$fq2_pref[$i]" . "$opts{'fq2_suff'}"; 
		&exeCmd_1cmd( "$opts{'exe_java'} -jar $opts{'jar_trimmo'} PE $opts{'paraJar_trimmo'}   $fq1   $fq2   $fq1_pref[$i]_pTr_R1.fq $fq1_pref[$i]_sTr_R1.fq $fq2_pref[$i]_pTr_R2.fq $fq2_pref[$i]_sTr_R2.fq $opts{'paraPE_trimmo'}" ); 
	} elsif ( $fq1_pref[$i] ne '' ) {
		&exeCmd_1cmd( "$opts{'exe_java'} -jar $opts{'jar_trimmo'} SE $opts{'paraJar_trimmo'}   $fq1_pref[$i]$opts{'fq1_suff'}   $fq1_pref[$i]_sTr_R1.fq   $opts{'paraSE_trimmo'}" ); 
	} elsif ( $fq2_pref[$i] ne '' ) {
		&exeCmd_1cmd( "$opts{'exe_java'} -jar $opts{'jar_trimmo'} SE $opts{'paraJar_trimmo'}   $fq2_pref[$i]$opts{'fq2_suff'}   $fq2_pref[$i]_sTr_R2.fq   $opts{'paraSE_trimmo'}" ); 
	} else {
		&tsmsg("[Wrn] Skip null line [$i].\n"); 
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

