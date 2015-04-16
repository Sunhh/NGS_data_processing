package fileSunhh; 
# 2014-11-22

use strict; 
use warnings; 
use File::Which; 
# use IO::Zlib; # This package doesn't support bzip2 files. And adding IO::Compress:Bzip2 and IO::Uncompress::Bunzip2 is too many work currently. Not very useful since I am using linux usually. 
use File::Basename; 
use File::Copy::Recursive; 
use File::Spec::Functions qw( catfile path ); 
use LogInforSunhh; 
use Cwd; 
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

=head1 renameByPat( [$file1, $file2, ...], $old_pattern, $new_pattern ) 

Function   : Rename files by pattern. 

=cut
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


=head1 openFH( $filename, $open_type )

Required   : $filename

Function   : Open file and return file handle one at a time. 

Input      : 
 $filename  : filename ended with .gz or .bz2 will be treated in compressed format. 
 $open_type : Defined in %goodFileType. Could be (read/write) + (|BZ2|GZ) 

=cut 
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

=head1 _chkExist( $file_to_check_existence )

Check if given path exists ( or in $PATH environment );

Mainly copy from http://stackoverflow.com/questions/8243189/check-if-file-exists-even-on-path

Return     : ($first_found_path, [$found_path_1, $found_path_2, ...])
=cut
sub _chkExist ($) {
	my $inFile = shift; 
	( defined $inFile and $inFile ne '' ) or return (undef(), []); 
	my @path = File::Spec::Functions::path; 
	my @pathext = ( q{} ); # This is for windows' extention names. 
	if ( $^O eq 'MSWin32' ) {
		push(@pathext, map { lc($_) } split(/;/, $ENV{PATHEXT})); 
	}
	my @results; 
	-e $inFile and push(@results, $inFile); 
	for my $dir ( @path ) {
		for my $ext ( @pathext ) {
			my $f = File::Spec::Functions::catfile( $dir, "$inFile$ext" ); 
			-e $f and push(@results, $f); 
		}
	}
	@results > 0 or return(undef(), []); 
	return($results[0], \@results); 
}# sub _chkExist () 

=head1 _dircopy()

Invoke File::Copy::Recursive::dircopy
=cut
sub _dircopy {
	return File::Copy::Recursive::dircopy(@_); 
}# sub _dircopy() 


1; 
