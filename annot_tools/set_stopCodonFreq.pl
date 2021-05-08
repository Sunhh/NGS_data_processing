#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
!@ARGV and die "perl $0 etrain_rawGenes_findErr.std augustus/config/species/orgName/orgName_parameters.cfg\n"; 

my $f1 = shift; 
my $f2 = shift; 


my %info = &info_from_etrain($f1); 
&replace_value($f2, \%info, [qw/freq_tag freq_taa freq_tga/]); 


sub info_from_etrain {
	my ($fn) = @_; 
	my %back; 
	open F,'<',"$fn" or die; 
	while (<F>) {
		m!^Frequency of stop codons! or next; 
		my $l_tag = <F>; 
		$l_tag =~ m!^tag:\s*\d+\s*\(([0-9.]+)\)\s*$! or die "$l_tag\n"; 
		$back{'freq_tag'} = $1; 
		my $l_taa = <F>; 
		$l_taa =~ m!^taa:\s*\d+\s*\(([0-9.]+)\)\s*$! or die "$l_taa\n"; 
		$back{'freq_taa'} = $1; 
		my $l_tga = <F>; 
		$l_tga =~ m!^tga:\s*\d+\s*\(([0-9.]+)\)\s*$! or die "$l_tga\n"; 
		$back{'freq_tga'} = $1; 
		last; 
	}
	close F; 
	return(%back); 
}# info_from_etrain() 

sub replace_value {
	my ($fn, $info_hr, $tag_ar) = @_; 
	my $backup_fnN = 1; 
	my $backup_fn = "$fn.old.$backup_fnN"; 
	while (-e $backup_fn) {
		$backup_fnN ++; 
		$backup_fn = "$fn.old.$backup_fnN"; 
		$backup_fnN >= 1000 and die "Existing $backup_fn\n"; 
	}
	my %toChangePara; 
	for my $k1 (@$tag_ar) {
		defined $info_hr->{$k1} or do { &tsmsg("[Wrn] Failed to find [$k1] from information.\n"); next; }; 
		$toChangePara{$k1} = 1; 
	}
	scalar(keys %toChangePara) > 0 or do { &tsmsg("[Wrn] No need to change file [$fn]\n"); return; }; 
	my %cfg2hs; 
	$cfg2hs{'/Constant/amberprob'} = 'freq_tag'; 
	$cfg2hs{'/Constant/ochreprob'} = 'freq_taa'; 
	$cfg2hs{'/Constant/opalprob'}  = 'freq_tga'; 
	my %hs2cfg = map { ($cfg2hs{$_}, $_); } keys %cfg2hs; 
	for my $k1 (keys %toChangePara) {
		defined $hs2cfg{$k1} or &tsmsg("[Wrn] No definition for [$k1]\n"); 
	}
	# $cfg2hs{''} = ''; 
	open F,'<',"$fn" or die;
	open O,'>',"$backup_fn" or die; 
	while (<F>) {
		print O $_; 
	}
	close O; 
	close F; 
	open F1,'<',"$backup_fn" or die; 
	open O1,'>',"$fn" or die; 
	while (<F1>) {
		chomp; 
		if      (m!^\s*#! or m!^\s*$!) {
			print O1 "$_\n"; 
		} elsif (s!^(\S+)(\s+)(\S+)!!) {
			my ($key_txt, $space_txt, $val_txt) = ($1, $2, $3); 
			my $hrID = (defined $cfg2hs{$key_txt}) ? $cfg2hs{$key_txt} : $key_txt; 
			defined $toChangePara{$hrID} and $val_txt = $info_hr->{$hrID}; 
			print O1 join('', $key_txt, $space_txt, $val_txt, $_)."\n"; 
		} else {
			&stopErr("[Err] Failed to parse line : $_\n"); 
		}
	}
	close O1; 
	close F1; 
	
}

