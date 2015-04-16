#!/usr/bin/perl
use strict; 
use warnings; 
use Getopt::Long; 
use fileSunhh; 
use LogInforSunhh; 

my %opts; 
GetOptions(\%opts, 
	"help!", 
	"augDir:s", # /data/Sunhh/src/Annot/maker/maker/exe/augustus
	"copy_species:s", 
); 


sub usage {
	print <<HH;
################################################################################
# perl $0 
# -help 
#
# Basic settings: 
# -augDir          [/data/Sunhh/src/Annot/maker/maker/exe/augustus/]
#                   If not given, 
#                   first check from \$AUGUSTUS_CONFIG_PATH
#                   then check if directory '/data/Sunhh/src/Annot/maker/maker/exe/augustus' exist. 
#                   Will report warn if not found any. 
#
# -copy_species     [SpeciesName1,SpeciesName2] Copy SpeciesName1 to SpeciesName2. 
################################################################################
HH
	exit 1; 
}

$opts{'help'} and &usage(); 
%opts == 0 and &usage(); 

# Basic settings. 
my %infor; 
$infor{'augDir'} = &_getAugDir(); 

# Invoke functiosn. 
&copySpecies() if ( defined $opts{'copy_species'} ); 


# Sub functions. 
sub copySpecies {
	my @spec = split(/,/, $opts{'copy_species'}); 
	s!\s!!g foreach (@spec); 
	@spec == 2 or &stopErr("[Err] Please use two species only.\n"); 
	&copySpec1to2($spec[0], $spec[1]); 
	return undef(); 
}

# Supporting functions. 
=head1 _getAugDir() 

First check $opts{'augDir'}, then check $ENV{'AUGUSTUS_CONFIG_PATH'}, 
then check /data/Sunhh/src/Annot/maker/maker/exe/augustus/, 
fail if none work. 
=cut
sub _getAugDir {
	defined $opts{'augDir'} and &fileSunhh::_chkExist($opts{'augDir'}) and return $opts{'augDir'}; 
	if ( defined $ENV{'AUGUSTUS_CONFIG_PATH'} and $ENV{'AUGUSTUS_CONFIG_PATH'} ne '' and -d "$ENV{'AUGUSTUS_CONFIG_PATH'}/../") {
		return "$ENV{'AUGUSTUS_CONFIG_PATH'}/../"; 
	} elsif ( -d '/data/Sunhh/src/Annot/maker/maker/exe/augustus/' ) {
		return '/data/Sunhh/src/Annot/maker/maker/exe/augustus/'; 
	} else {
		&stopErr("[Err] Please provide augustus directory -augDir\n"); 
	}
}# sub _getAugDir()

=head1 copySpec1to2($specFrom, $specTo, 'augCfgPath'=>"$augDir/config/")

Return        : 1 if failed. 
=cut
sub copySpec1to2 {
	my $sp1 = shift; 
	my $sp2 = shift; 
	my %parm = @_; 
	unless ( defined $parm{'augCfgPath'} ) {
		$parm{'augCfgPath'} = "$infor{'augDir'}/config/"; 
	}
	-d "$parm{'augCfgPath'}/species/$sp2/" and &stopErr("[Err] $parm{'augCfgPath'}/species/$sp2/ already exists.\n"); 
	fileSunhh::_dircopy( "$parm{'augCfgPath'}/species/$sp1/", "$parm{'augCfgPath'}/species/$sp2/" ); 
	my $fpath1 = &_speciesCfgFiles($sp1); 
	my $fpath2 = &_speciesCfgFiles($sp2); 
	for my $tk (keys %{$fpath1->{'base'}}) {
		my $f1    = $parm{'augCfgPath'} . '/species/' . $sp1 . '/' . $fpath1->{'base'}{$tk}; 
		my $f2    = $parm{'augCfgPath'} . '/species/' . $sp2 . '/' . $fpath2->{'base'}{$tk}; 
		my $f2old = $parm{'augCfgPath'} . '/species/' . $sp2 . '/' . $fpath1->{'base'}{$tk}; 
		-f $f2old or do { &tsmsg("[Wrn] $f2old not exist.\n"); next; }; 
		open F,'<',"$f2old" or die; 
		open O,'>',"$f2" or die; 
		while (<F>) {
			s!(\s)$sp1!$1$sp2!o; 
			print O $_; 
		}
		close O; close F; 
		unlink($f2old); 
	}
	for my $tk (keys %{$fpath1->{'more'}}) {
		my $f1    = $parm{'augCfgPath'} . '/species/' . $sp1 . '/' . $fpath1->{'more'}{$tk}; 
		my $f2    = $parm{'augCfgPath'} . '/species/' . $sp2 . '/' . $fpath2->{'more'}{$tk}; 
		my $f2old = $parm{'augCfgPath'} . '/species/' . $sp2 . '/' . $fpath1->{'more'}{$tk}; 
		-f $f2old or do { &tsmsg("[Wrn] $f2old not exist.\n"); next; }; 
		open F,'<',"$f2old" or die; 
		open O,'>',"$f2" or die; 
		while (<F>) {
			s!(\s)$sp1!$1$sp2!o; 
			print O $_; 
		}
		close O; close F; 
		unlink($f2old); 
	}
	return undef(); 
}# sub copySpec1to2()

=head1 _speciesCfgFiles($speciesName)

Return         : \%file_paths
 {'base'}{'...'} basic files from generic species.
 {'more'}{'...'} more files including .pbl files. 
=cut
sub _speciesCfgFiles {
	my $spName = shift; 
	my %paths; 
	$paths{'base'}{'cfgfilename'}     = $spName . '_parameters.cfg'; 
	$paths{'base'}{'weightfilename'}  = $spName . '_weightmatrix.txt'; 
	$paths{'base'}{'metafilename'}    = $spName . '_metapars.cfg'; 
	$paths{'base'}{'metautrfilename'} = $spName . '_metapars.utr.cfg'; 
	$paths{'base'}{'ovlpfilename'}    = $spName . '_ovlp_len.pbl'; 
	$paths{'base'}{'transshadowfilename'} = $spName . '_trans_shadow_bacterium.pbl'; 
	$paths{'more'}{'exonpblfilename'} = $spName . '_exon_probs.pbl'; 
	$paths{'more'}{'utrpblfilename'}  = $spName . '_utr_probs.pbl'; 
	$paths{'more'}{'intronpblfilename'} = $spName . '_intron_probs.pbl'; 
	$paths{'more'}{'igenicpblfilename'} = $spName . '_igenic_probs.pbl'; 
	# Anything else? 
	return \%paths; 
}# sub _speciesCfgPath() 


