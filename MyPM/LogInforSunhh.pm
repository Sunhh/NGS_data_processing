package LogInforSunhh; 
# Change output of exeCmd; 

use strict; 
use warnings; 
use Exporter qw(import);

our @EXPORT = qw(tsmsg stopErr exeCmd);
our @EXPORT_OK = qw();


# Print (to STDERR) input text with time stamp. 
sub tsmsg {
	my $tt = scalar( localtime() ); 
	print STDERR join('', "[$tt]", @_); 
}#End tsmsg() 

# Print (to STDERR) input text with time stamp, and then exit. 
sub stopErr {
	&tsmsg(@_); 
	exit (1); 
}#End stopErr() 

# Execute command with system and time record. 
sub exeCmd {
	for my $cmd ( @_ ) {
		&tsmsg("[CMD] $cmd\n"); 
		if ( system($cmd) == 0 ) {
			&tsmsg("[CMD_done]\n"); 
		} else {
			if ( $? == -1 ) {
				&tsmsg("[CMD_err] Failed to execute: $!\n"); 
			} elsif ( $? & 127 ) {
				&tsmsg("[CMD_err] Child died with signal ", $? & 127, ", ", ($? & 128)? 'with' : 'without', "coredump\n"); 
			} else {
				&tsmsg("[CMD_err] Child exited with value ", $? >> 8, "\n"); 
			}
		}
	}
}#End exeCmd 

1; 
