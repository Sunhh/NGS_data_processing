#!/usr/bin/perl
# Edit 20180507 : Add 'GN' value, but can't find 'SV' value
#  The rest differences contain '(strain ...)', 
use strict;

-t and !@ARGV and die "perl $0 uniprot_sprot_plants.dat > uniprot_sprot_plants.fasta\n"; 

my %hash; 

while (<>) {
	m/^\s*$/ and next; 
	chomp; 
	if (m!^\s*//\s*$!) {
		&output_hash(); 
		%hash = (); 
	} elsif ( m!^\s+! ) {
		$hash{'prev_tag'} eq 'SQ' or die "Bad line: $_\n$hash{'prev_tag'}\n"; 
		if ($hash{'prev_txt'} =~ m!;\s*$!) {
			$hash{'SQ'} .= $_; 
		} elsif (defined $hash{'SQ'}) {
			$hash{'SQ'} .= $_; 
		}
	} elsif ( m!^SQ\s{3}! ) {
		$hash{'prev_tag'} = 'SQ'; 
	} elsif ( m!^DE\s{3}(?:RecName|SubName)! ) {
		m!^DE\s{3}(?:RecName|SubName):\sFull=(.+);$! or die "Bad line: $_\n"; 
		if (!defined $hash{'DEf'}) {
			&set_h(\%hash, 'DEf', $1 ) and die "Repeat: Def; $_\n"; 
		} else {
		}
		$hash{'prev_tag'} = 'DE'; 
	} elsif ( m!^DE\s{3}Flags: (\S+)! ) {
		my $v = $1; 
		if ( $_ =~ m!(Fragment\S*)! ) {
			$v = $1; 
			&set_h(\%hash, 'FLG', $v) and die "Repeat: FLG; $_\n"; 
		}
	} elsif ( m!^ID\s{3}(\S+)! ) {
		&set_h(\%hash, 'ID', $1) and die "Repeat: ID; $_\n"; 
		$hash{'prev_tag'} = 'ID'; 
	} elsif ( m!^AC\s{3}([^\s;]+);! ) {
		if ( $hash{'prev_tag'} eq 'AC' and defined $hash{'AC'} ) {
		} else {
			&set_h(\%hash, 'AC', $1) and die "Repeat: AC; $_\n"; 
			$hash{'prev_tag'} = 'AC' ; 
		}
	} elsif ( m!^OS\s{3}(.+)$! ) {
		if ($hash{'prev_tag'} eq 'OS') {
			$hash{'OS'} .= " $1"; 
		} else {
			&set_h(\%hash, 'OS', $1) and die "Repeat: OS; $_\n"; 
			$hash{'prev_tag'} = 'OS'; 
		}
	} elsif ( m!^PE\s{3}(\d+):! ) {
		&set_h(\%hash, 'PE', $1) and die "Repeat: PE; $_\n"; 
		$hash{'prev_tag'} = 'PE'; 
	} elsif ( s!^GN\s{3}!! ) {
		if ($hash{'prev_tag'} eq 'GN') {
			$hash{'GN'} .= " $_"; 
		} else {
			&set_h(\%hash, 'GN', $_) and die "Repeat: GN; $_\n"; 
			$hash{'prev_tag'} = 'GN'; 
		}
	} elsif ( m!^DT\s{3}\S+,\s+sequence version (\d+)\.! ) {
		&set_h(\%hash, 'SV', $1) and die "Repeat: SV; $_\n"; 
		$hash{'prev_tag'} = 'SV'; 
	} else {
	}
	$hash{'prev_txt'} = $_; 
}

sub set_h {
	my ($h, $key, $val) = @_; 
	defined $h->{$key} and return(1); 
	$h->{$key} = $val; 
	return(); 
}

sub output_hash {
	(keys %hash == 0) and return(); 
	defined $hash{'ID'} or die "No ID found [$hash{'AC'}]\n"; 
	defined $hash{'AC'} or die "No AC found for [$hash{'ID'}]\n"; 
	defined $hash{'DEf'} or die "No DE RecName [$hash{'ID'}]\n"; 
	defined $hash{'OS'} or die "No OS found for [$hash{'ID'}]\n"; 
	defined $hash{'PE'} or die "No PE found for [$hash{'ID'}]\n"; 
	defined $hash{'SQ'} or die "No SQ for [$hash{'SQ'}]\n"; 
	my $save_SQ = ''; 
	while ($hash{'OS'} =~ s!\(([^()]*\([^()]*\)[^()]*|[^()]*)\)\s*\.?$!!g) {
		my $v = $1; 
		$v =~ m!^strain! and $save_SQ = " ($v)"; 
		$v eq "" and die "$hash{'OS'}\n$save_SQ\n$1\n"; 
	}
	$hash{'OS'} =~ s!\s+$!!; 
	$hash{'OS'} = "$hash{'OS'}$save_SQ"; 
	$hash{'OS'} !~ m!\bsp\.$! and $hash{'OS'} =~ s!\.$!!; 
	my $txt_GN = ''; 
	if (defined $hash{'GN'}) {
		$hash{'GN'} =~ s!\s*\{ECO:[^{}]+\}!!; 
		$hash{'GN'} =~ s!^[^\s;]+?=(.+?)\s*;($|\s.*$)!$1!; 
		$hash{'GN'} =~ s!\,\s.*$!!; 
		$txt_GN = " GN=$hash{'GN'}"; 
	}
	my $txt_SV = (defined $hash{'SV'}) ? " SV=$hash{'SV'}" : ''; 
	$hash{'DEf'} =~ s!\s*\{ECO:[^{}]+\}!!; 
	defined $hash{'FLG'} and $hash{'FLG'} =~ s!;+$!!; 
	my $txt_FLG = (defined $hash{'FLG'}) ? " ($hash{'FLG'})" : ''; 
	print STDOUT ">sp|$hash{'AC'}|$hash{'ID'} $hash{'DEf'}${txt_FLG} OS=$hash{'OS'}${txt_GN} PE=$hash{'PE'}${txt_SV}\n"; 
	$hash{'SQ'} =~ s!\s!!g; 
	$hash{'SQ'} =~ s!(.{60})!$1\n!g; chomp($hash{'SQ'}); 
	print STDOUT "$hash{'SQ'}\n"; 
}# output_hash() 
