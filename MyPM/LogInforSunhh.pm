package LogInforSunhh; 
# Change output of exeCmd; 

use strict; 
use warnings; 
use IPC::Open3; 
use Symbol; 
use Exporter qw(import);

our @EXPORT = qw(tsmsg stopErr exeCmd exeCmd_1cmd);
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
sub exeCmd_1cmd {
	my $cmd = shift; 
	my $is_print = shift; 
	$is_print //= 0; 
	
	if ( $is_print != 0 ) {
		&tsmsg("[CMD_print]$cmd\n"); 
		return; 
	}
	&tsmsg("[CMD] $cmd\n"); 
	if ( system($cmd) == 0 ) {
		&tsmsg("[CMD_done]$cmd\n"); 
	} else {
		if ( $? == -1 ) {
			&tsmsg("[CMD_err] Failed to execute: $!\n"); 
		} elsif ( $? & 127 ) {
			&tsmsg("[CMD_err] Child died with signal ", $? & 127, ", ", ($? & 128)? 'with' : 'without', "coredump\n"); 
		} else {
			&tsmsg("[CMD_err] Child exited with value ", $? >> 8, "\n"); 
		}
	}
}#End exeCmd 

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

### This is a method to monitor process in perl. 
### Copied from maker perl script, and edited minor. 
sub run {
	my $command = shift; 
	( defined $command and $command ne '' ) or return; 
	&tsmsg("[CMD] $command\n"); 
	my ($CHLD_IN, $CHLD_OUT, $CHLD_ERR) = (gensym, gensym, gensym); 
	my $pid = open3( $CHLD_IN, $CHLD_OUT, $CHLD_ERR, $command ); 
	{
		local $/ = \1; 
		while ( my $line = <$CHLD_ERR> ) {
			print STDERR $line; 
		}
	}
	waitpid $pid, 0; 
	&stopErr("[Err] Cmd failed: $command\n") if $? != 0; 
	return; 
}# End sub run. 


1; 
