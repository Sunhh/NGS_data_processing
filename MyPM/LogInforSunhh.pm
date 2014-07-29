package LogInforSunhh; 

use strict; 
use warnings; 
use Exporter qw(import);

our @EXPORT = qw(tsmsg stopErr);
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


1; 