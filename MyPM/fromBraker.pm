package fromBraker; 
use fileSunhh; 
use Exporter qw(import);
our @EXPORT = qw();
our @EXPORT_OK = qw();
#### Mainly copy functions from braker for analysing augustus result. 


############################################################
##  Sub-routines.
#############################################################

=head1 accuracy_calculator( august_test.stdout )

Calculate the result of testing AUGUSTUS on genbank files in a single number

Return : ( $accuracy_counted_from_all_classes )

=cut
sub accuracy_calculator{
	my $aug_out=shift;
	my ($nu_sen, $nu_sp, $ex_sen, $ex_sp, $gen_sen, $gen_sp);
	my $ifh = &openFH($aug_out, '<'); 
	while(<$ifh>){
		if(/^nucleotide level\s*\|\s*(\S+)\s*\|\s*(\S+)/){
			$nu_sen=$1;
			$nu_sp=$2;
		}
		if(/^exon level\s*\|.*\|.*\|.*\|.*\|.*\|\s*(\S+)\s*\|\s*(\S+)/){
			$ex_sen=$1;
			$ex_sp=$2;
		}
		if(/^gene level\s*\|.*\|.*\|.*\|.*\|.*\|\s*(\S+)\s*\|\s*(\S+)/){
			$gen_sen=$1;
			$gen_sp=$2;
		}
	}
	close($ifh); 
	my $target=(3*$nu_sen+2*$nu_sp+4*$ex_sen+3*$ex_sp+2*$gen_sen+1*$gen_sp)/15;
	return $target;
}# accuracy_calculator() 

=head1 setParInConfig( $fn_AUGUSTUS_species_parameters.cfg, $tag_to_change, $value_to_set )

WARNING!!! This will directly change the input file! 

##########################################
# change a parameter in a config file    #
# assume the format                      #
# parName    value   # comment           #
##########################################

=cut
sub setParInConfig{
	my $configFileName = shift;
	my $parName = shift;
	my $value = shift;
	open(CFGFILE, "+<$configFileName") or die ("Could not read config file $configFileName\n");
	my @lines = <CFGFILE>;
	foreach my $line (@lines){
		$line =~ s/(\s*$parName +)(\S+?)(\s|\#|$)/$1$value$3/;
	}
	seek(CFGFILE, 0,0);
	print CFGFILE @lines or die ("Could not write $configFileName");
	truncate(CFGFILE, tell(CFGFILE));
	close(CFGFILE);
}# setParInConfig()


1;

