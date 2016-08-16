package fileSunhh; 
# 2014-11-22

use strict; 
use warnings; 
use File::Which; 
# use IO::Zlib; # This package doesn't support bzip2 files. And adding IO::Compress:Bzip2 and IO::Uncompress::Bunzip2 is too many work currently. Not very useful since I am using linux usually. 
use File::Basename; 
use File::Copy::Recursive; 
use File::Copy; 
use File::Spec::Functions qw( catfile path ); 
use File::Path; 
use LogInforSunhh; 
use mathSunhh; 
use Cwd; 
use Exporter qw(import); 

# our @EXPORT = qw(tsmsg stopErr exeCmd);
our @EXPORT = qw(openFH renameByPat isSkipLine splitL wantLine wantLineC);
our @EXPORT_OK = qw();

my %goodFileType = qw(
	<           read
	>           write
	>>          add
	read        read
	write       write
	add         add
	<gz         readGZ
	>gz         writeGZ
	<bz2        readBZ2
	>bz2        writeBZ2
	readGZ      readGZ
	writeGZ     writeGZ
	readBZ2     readBZ2
	writeBZ2    writeBZ2
	addGZ       addGZ
	addBZ2      addBZ2
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

=head1 write2file( $filename, $text, $open_type )

Required   : $filename $text

Function   : Open file and write $text into file according to open_type; 

Return     : 0

=cut 
sub write2file {
	my ($fn, $txt, $open_type) = @_; 
	$open_type //= 'write'; 
	my $fh = &openFH($fn, $open_type); 
	print {$fh} $txt; 
	close ($fh); 
	return 0; 
}# sub write2file ()

=head1 load_bn6File ( $filename )

Input      : .bn6 file ( blastn -outfmt 6 )

Return     : (\@bn6Info)
     [idx]{'arr'} = [split(/\t/, $_)]
     [idx]{'k2v'}{$key} = $value; 
       $key could be qw/qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand/
=cut
sub load_bn6File {
	my ($fn) = @_; 
	my $inFh = &openFH($fn, '<'); 
	my @bn6Info; 
	while (&wantLineC($inFh)) {
		my @ta = &splitL("\t", $_); 
		$ta[0] eq 'qseqid' and next; 
		push(@bn6Info, {}); 
		$bn6Info[-1]{'arr'} = \@ta; 
		$bn6Info[-1]{'k2v'}{'qseqid'}    = $ta[0]; 
		$bn6Info[-1]{'k2v'}{'sseqid'}    = $ta[1]; 
		$bn6Info[-1]{'k2v'}{'pident'}    = $ta[2]; 
		$bn6Info[-1]{'k2v'}{'length'}    = $ta[3]; 
		$bn6Info[-1]{'k2v'}{'mismatch'}  = $ta[4]; 
		$bn6Info[-1]{'k2v'}{'gapopen'}   = $ta[5]; 
		$bn6Info[-1]{'k2v'}{'qstart'}    = $ta[6]; 
		$bn6Info[-1]{'k2v'}{'qend'}      = $ta[7]; 
		$bn6Info[-1]{'k2v'}{'sstart'}    = $ta[8]; 
		$bn6Info[-1]{'k2v'}{'send'}      = $ta[9]; 
		$bn6Info[-1]{'k2v'}{'evalue'}    = $ta[10]; 
		$bn6Info[-1]{'k2v'}{'bitscore'}  = $ta[11]; 
		$bn6Info[-1]{'k2v'}{'qlen'}      = $ta[12] // ''; 
		$bn6Info[-1]{'k2v'}{'slen'}      = $ta[13] // ''; 
		$bn6Info[-1]{'k2v'}{'sstrand'}   = $ta[14] // ''; 
		$bn6Info[-1]{'k2v'}{'bitscore'} =~ s!^\s+|\s+$!!g; 
	}
	return (\@bn6Info); 
}# load_bn6File() 



=head1 load_tabFile ( $filename , $keep_annot<0> ) 

Function   : Load in a tab delimited file into a array, whose element is a array reference of all files of the line. 
             If $keep_annot is TRUE (1), those lines beginning with '#' or blank lines will also be loaded in. 

Return     : ( $line_1_array, $line_2_array, ... )

  $line_1_array = [ field_1, field_2, ... ]

=cut
sub load_tabFile {
	my ($fn, $keep_annot) = @_; 
	$keep_annot //= 0; 
	my $fh = &openFH($fn, '<'); 
	my @back; 
	if ( $keep_annot ) {
		while (<$fh>) {
			$_ =~ s/[\r\n]+$//; 
			push(@back, [ &splitL("\t", $_) ]); 
		}
	} else {
		while (&wantLineC($fh)) {
			push(@back, [ &splitL("\t", $_) ]); 
		}
	}
	close($fh); 
	return(@back); 
}# load_tabFile ()



=head1 load_agpFile( $filename , $if_save_unknown_str [0|1])

Required   : $filename

Function   : load in an agp file into a hash, which can be input into $mathSunhh_obj->switch_position(); 

Return     : ( \%ctg2scf )

  %ctg2scf = ( $ctgID => [ [ctgS, ctgE, scfID, scfS, scfE, scfStr(+/-/?)], [], ... ] ) # This is sorted. 
  If the input $if_save_unknown_str is 0, the scfStr will be either + or -. 

=cut
sub load_agpFile {
	my ( $fn , $save_str ) = @_; 
	$save_str //= 0; 
	my $fh = &openFH( $fn, '<' ); 
	my %ctg2scf; 
	while (<$fh>) {
		m/^\s*(#|$)/ and next; 
		chomp; s/[^\S\t]+$//; 
		my @ta = &splitL("\t", $_); 
		$ta[4] eq 'W' or next; 
		$save_str or ( $ta[8] eq '?' and $ta[8] = '+' ); 
		push( @{$ctg2scf{$ta[5]}}, [@ta[ 6,7,0,1,2,8 ]] ); 
	}
	close($fh); 
	for my $tk (keys %ctg2scf) {
		@{$ctg2scf{$tk}} = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{ $ctg2scf{$tk} }; 
	}
	
	return(\%ctg2scf); 
}# load_agpFile ()


=head1 openFH( $filename, $open_type )

Required   : $filename

Function   : Open file and return file handle one at a time. 
 $open_type could be '<|>|read|write|add|readGZ|writeGZ|readBZ2|writeBZ2|addGZ|addBZ2'

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
	} elsif ($type eq 'addGZ' or $type eq 'addBZ2') { 
		$fh = &oCompressFile($f, $type); 
	} elsif ($type eq 'add') {
		$f =~ m/\.(gz|gzip)$/   and do { $type = "addGZ"; &tsmsg("[Msg] Using gzip format for [$f]\n"); goto RE_CHK; }; 
		$f =~ m/\.(bz2|bzip2)$/ and do { $type = "addBZ2"; &tsmsg("[Msg] Using bzip2 format for [$f]\n"); goto RE_CHK; }; 
		open ($fh, '>>', "$f") or &stopErr("[Err] $! [$f]\n"); 
	} else {
		# Something is wrong. 
		&stopErr("[Err]Unknown type[$type].\n"); 
	}
	return $fh; 
}#End openFH() 

=head1 splitL( "\t", $line )

Return         : (@splitted_line)

=cut
sub splitL {
	chomp($_[1]); 
	my @back = split(/$_[0]/, "$_[1]\n"); 
	chomp($back[-1]); 
	return (@back); 
}# splitL() 

=head1 wantLine( $fileHandle )

Return         : ($line_text)

Description    : This $line_text is not modified, but this line should not pass check of &isSkipLine(). 

=cut
sub wantLine {
	if (defined $_[0]) {
		while ( readline($_[0]) ) {
			&isSkipLine($_) and next; 
			return($_); 
		}
	} else {
		while ( readline() ) {
			&isSkipLine($_) and next; 
			return($_); 
		}
	}
	return; 
}# wantLine() 

=head1 wantLineC( $fileHandle )

Return         : ($line_text_wo_return)

Function       : Remove the tailing '\r|\n's in good input line and return it. 

=cut
sub wantLineC {
	$_ = &wantLine( @_ ); # Here and elsewhere I continuously use '$_' to ignore problem in using. '$_' is always changed when using, so do be careful with it. 
	defined $_ or return(); 
	$_ =~ s/[\r\n]+$//; 
	return( $_ ); 
}# wantLineC() 


=head1 log_section( $cntN_to_chk, \%tmp_cnt )

Return         : (1|0) # 1 - get to a new section divided by 'cntN_step|base' ; 0 - Not yet. 

Function       : %tmp_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>5e6 ); 
                 Edit %tmp_cnt if needed when revoked. 
                 This function is used for recording usage. 
                 When 'cntN_step' <= 0, always return 0. 
                 When $cntN_to_chk > $tmp_cnt{'cntN_base'} , return 1 and increase 'cntN_base'; 
=cut
sub log_section {
	my ($v, $hr) = @_; 
	$hr->{'cntN_step'} //= 5e6 ; 
	$hr->{'cntN_base'} //= 0 ; 
	$hr->{'cntN_step'} > 0 or return 0; 
	# $hr->{'cntN_step'} > 0 or &stopErr("[Err] 'cntN_step' for log_section() should > 0!\n"); 
	$v > $hr->{'cntN_base'} or return 0; 
	while ($v > $hr->{'cntN_base'} ) {
		$hr->{'cntN_base'} += $hr->{'cntN_step'}; 
	}
	return 1; 
}# log_section() 


=head1 isSkipLine( $line )

Return         : 0/1
   0 - for good lines 
   1 - for bad line which is empty or begins with '#' 

=cut
sub isSkipLine {
	$_[0] =~ m/^\s*(?:#|$)/ and return 1; 
	return 0; 
}# isSkipLine() 

=head1 new_tmp_dir( 'maxNum' => 999999 )

Function       : Creat a new temporary directory with name 'tmp\d+', and return the name of dir. 
 I will test from tmp0 to tmp$maxNum until the filename does not exist. 
 If all files exists, I will return undef() and report an warning. 
=cut
sub new_tmp_dir {
	my %parm = @_; 
	$parm{'maxNum'} //= 999999; 
	for (my $i=0; $i<=$parm{'maxNum'}; $i++) {
		my $fname = "tmp$i"; 
		-e $fname and next; 
		return $fname; 
	}
	&tsmsg("[Wrn] All tmp dir from tmp0 to tmp$parm{'maxNum'} exists!\n"); 
	return undef(); 
}# sub new_tmp_dir() 

=head1 new_tmp_file( 'maxNum' => 999999 )

Function       : Creat a new temporary file with name 'tmp\d+', and return the name of file. 
 I will test from tmp0 to tmp$maxNum until the filename does not exist. 
 If all files exists, I will return undef() and report an warning. 
=cut
sub new_tmp_file {
	my %parm = @_; 
	$parm{'maxNum'} //= 999999; 
	for (my $i=0; $i<=$parm{'maxNum'}; $i++) {
		my $fname = "tmp$i"; 
		-e $fname and next; 
		return $fname; 
	}
	&tsmsg("[Wrn] All tmp file from tmp0 to tmp$parm{'maxNum'} exists!\n"); 
	return undef(); 
}# sub new_tmp_file() 

=head1 dvd_file( $in_filename_or_filehandle, $num_of_groups, 'keep_order'=>0, 'with_header'=>0, 'sub_pref'=>'sub_' )

Function       : divide $in_filename_or_filehandle into subgroup_files within format 'sub_pref'$grpNumber . 
  If the line number of in_file is smaller than $num_of_groups, $num_of_in_file_lines (-1 if with header) sub_files will be created. 

Return         : ($sub_filename_1, $sub_filename_2, ...)
  If in_file is : 
header
1
2
3
4
  &dvd_file( $in_file, 2, 'keep_order'=>1, 'with_header'=>1, 'sub_pref'=>'out_' ) will return ('out_0', 'out_1'), and in the files: 
# out_0 contains: 
header
1
2
# out_1 contains: 
header
3
4

=cut
sub dvd_file {
	my $inFn = shift; 
	my $grpN = shift; 
	my $inFh; 
	if ( ref($inFn) eq '' ) {
		$inFh = &openFH( $inFn, '<' ); 
	} elsif ( ref($inFn) eq 'GLOB' ) {
		$inFh = $inFn; 
	} else {
		&stopErr("[Err] The 1st input [$inFn] of dvd_file() should be a filename or file_handle.\n"); 
	}
	( defined $grpN and int($grpN) > 0 ) or &stopErr("[Err] 2nd input [$grpN] of dvd_file() must not be smaller than 1.\n"); 
	$grpN = int($grpN); 
	my %parm = &mathSunhh::_setHashFromArr(@_); 
	$parm{'keep_order'} //= 0; 
	$parm{'with_header'} //= 0; 
	$parm{'sub_pref'} //= 'sub_'; 
	my $header = ''; 
	my @in_data = <$inFh>; 
	$parm{'with_header'} and $header = shift(@in_data); 
	my @sub_data = @{ &mathSunhh::dvd_array(\@in_data, $grpN, $parm{'keep_order'}) }; 
	
	my @sub_file_names; 
	for (my $i=0; $i<@sub_data; $i++) {
		my $subFn = "$parm{'sub_pref'}$i"; 
		push(@sub_file_names, $subFn); 
		&write2file( $subFn, join('',$header,@{$sub_data[$i]}), '>' ); 
	}

	return @sub_file_names; 
}# dvd_file() 



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
	} elsif ( $type eq 'addGZ' ) {
		$has_gzip ne '' or &stopErr("[Err] No gzip command found in PATH.\n"); 
		open($fh,'|-',"$has_gzip >> $f") or &stopErr("[Err] $! [$f]\n"); 
	} elsif ( $type eq 'addBZ2' ) {
		$has_bzip2 ne '' or &stopErr("[Err] No bzip2 command found in PATH.\n"); 
		open($fh,'|-',"$has_bzip2 >> $f") or &stopErr("[Err] $! [$f]\n"); 
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
=head1 _abs_path()

Invoke Cwd::abs_path() 
=cut
sub _abs_path {
	return Cwd::abs_path(@_); 
}# sub _abs_path() 
=head1 _copy()
Invoke File::Copy::copy()
=cut
sub _copy {
	return File::Copy::copy(@_); 
}#sub _copy()
=head1 _move()
Invoke File::Copy::move()
=cut
sub _move {
	return File::Copy::move(@_); 
}#sub _move()

=head1 _rmtree()
=cut
sub _rmtree {
	return File::Path::rmtree(@_); 
}#sub _rmtree()

=head1 _basename()
=cut
sub _basename {
	return File::Basename::basename(@_); 
}
=head1 _dirname()
=cut
sub _dirname {
	return File::Basename::dirname(@_); 
}
1; 
