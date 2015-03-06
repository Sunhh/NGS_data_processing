package fileSunhh; 
# 2014-11-22

use strict; 
use warnings; 
use File::Which; 
# use IO::Zlib; # This package doesn't support bzip2 files. And adding IO::Compress:Bzip2 and IO::Uncompress::Bunzip2 is too many work currently. Not very useful since I am using linux usually. 
use LogInforSunhh; 
use Exporter qw(import);

# our @EXPORT = qw(tsmsg stopErr exeCmd);
our @EXPORT = qw(openFH renameByPat);
our @EXPORT_OK = qw();

my %goodFileType = qw(
	<           read
	>           write
	read        read
	write       write
	<gz         readGZ
	>gz         writeGZ
	<bz2        readBZ2
	>bz2        writeBZ2
	readGZ      readGZ
	writeGZ     writeGZ
	readBZ2     readBZ2
	writeBZ2    writeBZ2
); 

# Check if there is gzip/bzip2 software in current system. 
my ($has_gzip, $has_bzip2) = (File::Which::which("gzip"), File::Which::which("bzip2")); 

# Rename files by pattern. 
# &renameByPat( [$file1, $file2, ...], $old_pattern, $new_pattern )
#  
sub renameByPat {
	my ($infile, $pat_old, $pat_new) = @_; 
	my @inFiles; 
	if ( ! (ref($infile)) ) {
		push(@inFiles, $infile); 
	} elsif ( ref($infile) eq 'ARRAY' ) {
		push(@inFiles, @$infile); 
	} elsif ( ref($infile) eq 'SCALAR' ) {
		push(@inFiles, $$infile); 
	} else {
		&stopErr("[Err] Input file [$infile] not recognized!\n"); 
	}
	
	my @backList; # ( [old_name_1, new_name_1], [old_name_2, new_name_2], ... )
	for my $infileID ( @inFiles ) {
		defined $infileID and -e $infileID  or do { &tsmsg("[Err] input file [$infileID] not found. Skipping\n"); next; }; 
		my $outfileID = $infileID; 
		$outfileID =~ s!$pat_old!$pat_new!; 
		&tsmsg("[Msg] Rename file from [$infileID] to [$outfileID] .\n"); 
		rename($infileID, $outfileID); 
		push(@backList, [$infileID, $outfileID]); 
	}
	return \@backList; 
}# sub renameByPat() 


# Open file and return file handle one at a time. 
# &openFH( $file_to_handle, $handle_type)
#   "file_to_handle" : filename ended with .gz or .bz2 will be treated in compressed format. 
#   "handle_type"    : Defined in %goodFileType. Could be (read/write) + (|BZ2|GZ)
sub openFH ($$) {
	my $f = shift; 
	my $type = shift; 
	defined $f or &stopErr("[Err] file name not defined in openFH()\n"); 
	defined $type or $type = 'read'; 
	defined $goodFileType{$type} or &stopErr("[Err]Unknown open method tag [$type].\n"); 
	$type = $goodFileType{$type}; 
	my $fh; 
	RE_CHK: 
	if ($type eq 'read') {
		$f =~ m/\.(gz|gzip)$/ and do { $type = "readGZ"; &tsmsg("[Msg] Using gzip format for [$f]\n"); goto RE_CHK; }; 
		$f =~ m/\.(bz2|bzip2)$/ and do { $type = "readBZ2"; &tsmsg("[Msg] Using bzip2 format for [$f]\n"); goto RE_CHK; }; 
		open ($fh, '<', "$f") or &stopErr("[Err] $! [$f]\n"); 
	} elsif ($type eq 'write') {
		$f =~ m/\.(gz|gzip)$/ and do { $type = "writeGZ"; &tsmsg("[Msg] Using gzip format for [$f]\n"); goto RE_CHK; }; 
		$f =~ m/\.(bz2|bzip2)$/ and do { $type = "writeBZ2"; &tsmsg("[Msg] Using bzip2 format for [$f]\n"); goto RE_CHK; }; 
		open ($fh, '>', "$f") or &stopErr("[Err] $! [$f]\n"); 
	} elsif ($type eq 'readGZ') {
		$fh = &iCompressFile($f, 'gz'); 
	} elsif ($type eq 'readBZ2') {
		$fh = &iCompressFile($f, 'bz2'); 
	} elsif ($type eq 'writeGZ') {
		$fh = &oCompressFile($f, 'gz'); 
	} elsif ($type eq 'writeBZ2') {
		$fh = &oCompressFile($f, 'bz2'); 
	} else {
		# Something is wrong. 
		&stopErr("[Err]Unknown type[$type].\n"); 
	}
	return $fh; 
}#End openFH() 

sub iCompressFile ($$) {
	my $f = shift; 
	my $type = shift; 
	my $fh; 
	if ( $type eq 'gz' ) {
		$has_gzip ne '' or &stopErr("[Err] No gzip command found in PATH.\n"); 
		open($fh,'-|',"$has_gzip -cd $f") or &stopErr("[Err] $! [$f]\n"); 
	} elsif ( $type eq 'bz2' ) {
		$has_bzip2 ne '' or &stopErr("[Err] No bzip2 command found in PATH.\n"); 
		open($fh,'-|',"$has_bzip2 -cd $f") or &stopErr("[Err] $! [$f]\n"); 
	} else {
		&stopErr("[Err] Unknown type [$type]\n"); 
	}
	return $fh; 
}#End iCompressFile()

sub oCompressFile ($$) {
	my $f = shift; 
	my $type = shift; 
	my $fh; 
	if ( $type eq 'gz' ) {
		$has_gzip ne '' or &stopErr("[Err] No gzip command found in PATH.\n"); 
		open($fh,'|-',"$has_gzip > $f") or &stopErr("[Err] $! [$f]\n"); 
	} elsif ( $type eq 'bz2' ) {
		$has_bzip2 ne '' or &stopErr("[Err] No bzip2 command found in PATH.\n"); 
		open($fh,'|-',"$has_bzip2 > $f") or &stopErr("[Err] $! [$f]\n"); 
	} else {
		&stopErr("[Err] Unknown type [$type]\n"); 
	}
	return $fh; 
}#End oCompressFile()

1; 
