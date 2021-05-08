#!/usr/bin/perl
# 20190717 : Replace database file in .asn.1 files to use blast_formatter. 
use strict; 
use warnings; 

!@ARGV and die "perl $0 new_db_path input_raw.asn > new.asn\n"; 

my $new_db = shift; 

while (<>) {
	if ( m!^\s+subject database \"! ) {
		s!^(\s+subject database \").+\"!$1${new_db}"!o or die "Failed at line: $_\n"; 
	}
	print ; 
}
