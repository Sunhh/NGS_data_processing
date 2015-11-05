#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use fastaSunhh; 
my $fas_obj = fastaSunhh->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"agpFile:s", 
	"Nlen:i", # 60
	"Ulen:i", # 60
); 
$opts{'Nlen'} //= 60; 
$opts{'Ulen'} //= 60; 

my $help_txt = <<HH; 

perl $0 -agpFile SP_Nov_11_13_30_8_90_3_superscaffold_merged_edit.agp 

-help 

-Nlen               [$opts{'Nlen'}] Length for 'N'
-Ulen               [$opts{'Ulen'}] Length for 'U'

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
for my $tk (qw/agpFile/) {
	defined $opts{$tk} or &LogInforSunhh::usage($help_txt); 
}

&tsmsg("[Msg] Loading AGP file [$opts{'agpFile'}].\n"); 
my @agp_arr = @{ &load_agp_file( $opts{'agpFile'} ) }; 

&tsmsg("[Msg] Output agp file.\n"); 
for my $la (@agp_arr) {
	if ( $la->[4] =~ m/^N$/i ) {
		$la->[5] = $opts{'Nlen'}; 
	} elsif ( $la->[4] =~ m/^U$/i ) {
		$la->[5] = $opts{'Ulen'}; 
	} elsif ( $la->[4] eq 'W' ) {
		; 
	} else {
		&stopErr("[Err] No rule for component_type[5] [$la->[4]] yet.\n"); 
	}
	print STDOUT join("\t", @$la)."\n"; 
}

&tsmsg("[Msg] Done $0.\n"); 

# [Sunhh@whale super_scaffold]$ head -3 SP_Nov_11_13_30_8_90_3_superscaffold_merged.agp
# ##agp-version   2.0
# Super_scaffold_248      1       341379  1       W       SpoScf_00639    1       341379  -
# Super_scaffold_248      341380  600696  2       N       259317  scaffold        yes     map
sub load_agp_file {
	# https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/ 
	my $fh = &openFH($_[0], '<'); 
	my @back; 
	while (<$fh>) {
		m/^\s*(#|$)/ and next; 
		chomp; 
		my @ta = split(/\t/, $_); 
		push(@back, [@ta]); 
	}
	close($fh); 
	return(\@back); 
}

