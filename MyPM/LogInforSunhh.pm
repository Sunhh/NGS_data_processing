package LogInforSunhh; 
# Change output of exeCmd; 

use strict; 
use warnings; 
use IPC::Open3; 
use Symbol; 
use Exporter qw(import);

our @EXPORT = qw(tsmsg stopErr exeCmd exeCmd_1cmd);
our @EXPORT_OK = qw();

=head1 tsmsg($seg1, $seg2, ...)

Function : All segN will be joined with empty character. 
           Print (to STDERR) input text with time stamp. 

=cut
sub tsmsg {
	my $tt = scalar( localtime() ); 
	print STDERR join('', "[$tt]", @_); 
}#End tsmsg() 

=head1 stopErr()

Function : Print (to STDERR) input text with time stamp, and then exit. 

=cut
sub stopErr {
	&tsmsg(@_); 
	exit (1); 
}#End stopErr() 

=head1 exeCmd_1cmd($cmd_to_exec, $print_only)

Function : Execute command with system and time record. 
           When $print_only == 1, will only print command instead of executing it. 

=cut
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
		return 0; 
	} else {
		if ( $? == -1 ) {
			&tsmsg("[CMD_err] Failed to execute: $! . CMD: $cmd\n"); 
		} elsif ( $? & 127 ) {
			&tsmsg("[CMD_err] Child died with signal ", $? & 127, ", ", ($? & 128)? 'with' : 'without', " coredump. CMD: $cmd\n"); 
		} else {
			&tsmsg("[CMD_err] Child exited with value ", $? >> 8, " . CMD: $cmd\n"); 
		}
	}
	return 1; 
}#End exeCmd 

=head1 timed_out( $cmd, $num_secs_to_timeout )

http://stackoverflow.com/questions/1962985/how-can-i-timeout-a-forked-process-that-might-hang
http://www.ahlinux.com/perl/4677.html

=cut
sub timed_out {
	my $retval ; 
	my ($cmd, $num_secs_to_timeout) = @_; 
	my $pid = fork; 
	if ( $pid > 0 ) { # This is parent process 
		eval {
			local $SIG{'ALRM'} = 
			  sub { kill 9, -$pid; &tsmsg("[CMD_timeout_err]$cmd\n"); $retval = 124; }; 
			alarm $num_secs_to_timeout; 
			waitpid( $pid, 0 ); 
			alarm 0; 
		}; 
		return defined($retval) ? $retval : 0 ; 
	} elsif ( $pid==0 ) { # This is child process 
		setpgrp(0,0); 
		&exeCmd_1cmd($cmd); 
	} else { # The fork is not successful. 
		&stopErr("Failed to exec fork.\n"); 
	}
	return; 
}# sub timed_out()


=head1 exeCmd($cmd1, $cmd2, $cmd3, ... )

Function : Execute commands with system and time record. 

=cut 
sub exeCmd {
	for my $cmd ( @_ ) {
		&tsmsg("[CMD] $cmd\n"); 
		if ( system($cmd) == 0 ) {
			&tsmsg("[CMD_done]$cmd\n"); 
			return 0; 
		} else {
			if ( $? == -1 ) {
				&tsmsg("[CMD_err] Failed to execute: $! . CMD: $cmd\n"); 
			} elsif ( $? & 127 ) {
				&tsmsg("[CMD_err] Child died with signal ", $? & 127, ", ", ($? & 128)? 'with' : 'without', " coredump. CMD: $cmd\n"); 
			} else {
				&tsmsg("[CMD_err] Child exited with value ", $? >> 8, " . CMD: $cmd\n"); 
			}
		}
	}
	return 1; 
}#End exeCmd 

=head1 LogInforSunhh::run($cmd)

### This is a method to monitor process in perl. 
### Copied from maker perl script, and edited minor. 

=cut
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

=head1 LogInforSunhh::usage($help_txt)

Print the $help_txt and then exit(1); 

=cut
sub usage {
	print STDERR $_[0]; 
	exit(1); 
}

1; 
