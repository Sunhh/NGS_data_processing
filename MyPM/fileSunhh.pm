package fileSunhh; 

use strict; 
use warnings; 
use IO::Zlib; 
use LogInforSunhh; 
use Exporter qw(import);

# our @EXPORT = qw(tsmsg stopErr exeCmd);
our @EXPORT = qw(openFH);
our @EXPORT_OK = qw();

my %goodFileType = qw(
	<	read
	>	write
	read	read
	write	write
	<gz	readGZ
	>gz	writeGZ
	<bz	readBZ
	>bz	writeBZ
	readGZ	readGZ
	writeGZ	writeGZ
	readBZ	readBZ
	writeBZ	writeBZ
); 

# Check if there is gzip/bzip2 software in current system. 
my ($has_gzip, $has_bzip2) = (0,0); 


# Open file and return file handle one at a time. 
sub openFH ($$) {
	my $f = shift; 
	my $type = shift; 
	defined $type or $type = 'read'; 
	defined $goodFileType{$type} or &stopErr("[Err]Unknown open method tag [$type].\n"); 
	local *FH; 
	if ($type eq 'read') {
	} elsif ($type eq 'write') {
	} elsif ($type eq 'readGZ') {
	} elsif ($type eq 'readBZ') {
	} elsif ($type eq 'writeGZ') {
	} elsif ($type eq 'writeBZ') {
	} else {
		# Something is wrong. 
		&stopErr("[Err]Unknown type[$type].\n"); 
	}
	return *FH; 
}#End openFH() 


1; 
