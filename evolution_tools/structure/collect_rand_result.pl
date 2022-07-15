#!/usr/bin/perl
use strict; 
use warnings; 
use File::Copy; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"rmLoci!", # Remove Locus lines; 
	"help!", 
); 

!@ARGV and die "perl $0 in_dir out_dir\n"; 


my $dir1 = shift; 
my $tgtDir = shift; 
$tgtDir //= 'collected_struct_result_f'; 
unless ( -e $tgtDir ) {
	mkdir($tgtDir); 
}
my $dir_base = $dir1; 
$dir_base =~ s!\/+$!!; 
$dir_base =~ s!^.*\/(\S+?)$!$1!; 

my %h; 

opendir(D1, $dir1); 
my @dd1 = grep { m!^structure_K(\d+)$! } readdir(D1); 
for my $td1 (@dd1) {
	if ( -e "$dir1/$td1/struct.result_f" ) {
		if ($opts{'rmLoci'}) {
			&filter_file( "$dir1/$td1/struct.result_f", "$tgtDir/${dir_base}.$td1.struct.result_f" ); 
		} else {
			File::Copy::copy("$dir1/$td1/struct.result_f", "$tgtDir/${dir_base}.$td1.struct.result_f"); 
		}
	} else {
		warn "No result found in $dir1/$td1/\n"; 
	}
}
closedir (D1); 

sub filter_file {
	my ($fi, $fo) = @_; 
	open FI,'<',"$fi" or die "$fi\n"; 
	open FO,'>',"$fo" or die "$fo\n"; 
	my %th1; 
	while (<FI>) {
		if (m!^Locus\s+\d+!) {
			$th1{'is_loci'} = 1; 
		} elsif (m!^\s*$!) {
			if ($th1{'is_loci'}) {
				%th1 = (); 
				next; 
			}
			%th1 = (); 
		}
		defined $th1{'is_loci'} and $th1{'is_loci'} == 1 and next; 
		print FO; 
	}
	close FO; 
	close FI; 
}

