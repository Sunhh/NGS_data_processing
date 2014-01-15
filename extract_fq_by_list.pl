#!/usr/bin/perl -w 
# 2013-11-13 Made for pick fastq read from list. 
use strict; 
use warnings; 
use Getopt::Long; 

my %opts; 

GetOptions(\%opts, 
	"help!", 
	"refLis:s", 
	"refFq:s", 
	"srcFq:s", 
	"outFq:s", 
	"rdKey!", 
	"trim12!", 
	"mode:s", # keep/drop
); 

sub usage {
	print STDOUT <<INFO; 
perl $0 -refLis in_list1,in_list2 -srcFq source1.fq,source2.fq -outFq source.fq.filt
-mode    keep/drop
-refFq   in_fq1,in_fq2
-rdKey   [!]
-trim12  [remove R1R2 tag]
INFO
	exit 1; 
}
if ($opts{help}) {
	&usage(); 
} elsif ( !defined $opts{refLis} and !defined $opts{refFq} ) {
	&usage(); 
# } elsif ( !defined $opts{} ) {
}


sub tmsg {
	my $txt = join('', @_); 
	my $tt = scalar( localtime() ); 
	print STDERR "[$tt]$txt"; 
}

my %has_tag; 
if (defined $opts{refLis}) {
	for my $fname ( split(/,/, $opts{refLis}) ) { 
		my $ls_fh = &openFH( $fname ); 
		&tmsg("[Rec] Reading refLis [$fname]\n"); 
		my $num_ls_lines=0;
		while (<$ls_fh>) {
			chomp; 
			my ($tk) = (split(/\t/, $_))[0]; 
			$has_tag{$tk} = 1; 
			if ( $opts{rdKey} ) {
				my ($tid = $tk) =~ s!^(\S+)\s.*!$1!; 
				$has_tag{$tid} = 1; 
			}
			if ( $opts{trim12} ) {
				my ($tid = $tk) =~ s!^(\S+)/[12](?:\s.*|)$!$1!; 
				$has_tag{$tid} = 1; 
			}
			$num_ls_lines ++; 
		}
		my $numInLis = scalar( keys %has_tag ); 
		close($ls_fh); 
		&tmsg( "[Rec] There are $num_ls_lines lines and $numInLis tags.\n" ); 
	}#End for 
}

if (defined $opts{refFq}) {
	for my $fname ( split(/,/, $opts{refFq}) ) { 
		my $fq_fh = &openFH( $fname ); 
		&tmsg("[Rec] Reading refFq [$fname]\n"); 
		my $num_rd_total=0; 
		my $exe_code_refFq = <<'FQ'; 
		while ( !eof($fq_fh) ) {
			$num_rd_total % 1e6 == 1 and &tmsg( "[Msg] $num_rd_total reads tagged.\n" ); 
			my ($fq_recR, $fq_idx, $fq_rdID) = &readRecord($fq_fh); 
			$has_tag{$fq_rdID} = 1; 
			$num_rd_total ++; 

FQ
		if ( $opts{rdKey} ) {
			$exe_code_refFq .= <<'FQ'; 
			( my $fq_rdID_key = $fq_rdID ) =~ s/^(\S+)\s.*/$1/; 
			$has_tag{$fq_rdID_key} = 1; 
FQ
		}
		if ( $opts{trim12} ) {
			$exe_code_refFq .= <<'FQ'; 
			( my $fq_rdID_trim = $fq_rdID ) =~ s!^(\S+)/[12](?:\s.*|)$!$1!; 
			$has_tag{$fq_rdID_trim} = 1; 
FQ
		}
		$exe_code_refFq .= <<'FQ'; 
		}# End while 
FQ
		eval("$exe_code_refFq"); 
		close ($fq_fh); 
		my $numInLis = scalar( keys %has_tag ); 
		&tmsg( "[Rec] There are $num_rd_total reads and $numInLis tags.\n" ); 
	}# End for 
}#End if refFq

defined $opts{outFq} or $opts{outFq} = ''; 
my ($tn, $num_keep, $num_drop) = (0, 0, 0); 
&tmsg( "[Rec] Out fastq file is [$opts{outFq}]\n" ); 
my $out_fh = ( $opts{outFq} eq '' )
	? \*STDOUT 
	: &openFH( $opts{outFq}, '>' )
; 

for my $fname ( split(/,/, $opts{srcFq}) ) { 
	&tmsg( "[Rec] Reading srcFq [$fname]\n" ); 
	my $src_fh = &openFH( $fname ); 

	my $exe_code_src = <<'SRC'; 
	while ( !eof($src_fh) ) {
		$tn % 1e6 == 1 and &tmsg("[Msg] $tn reads processed.\n"); 
		my ($fq_recR, $fq_idx, $fq_rdID) = &readRecord($src_fh); 

SRC

	my 	@dd_comp = ('defined $has_tag{$fq_rdID}'); 
	if ( $opts{rdKey} ) {
		$exe_code_src .= <<'SRC'; 
		( my $fq_rdID_key = $fq_rdID ) =~ s/^(\S+)\s.*/$1/; 
SRC
		push(@dd_comp, 'defined $has_tag{$fq_rdID_key}'); 
	}
	if ( $opts{trim12} ) {
		$exe_code_src .= <<'SRC'; 
		( my $fq_rdID_trim = $fq_rdID ) =~ s!^(\S+)/[12](?:\s.*|)$!$1!; 
SRC
		push(@dd_comp, 'defined $has_tag{$fq_rdID_trim}'); 
	}

	my $DD_comp = join(" or ", @dd_comp); 
	my $slct = (defined $opts{mode}) 
		? ( $opts{mode} eq 'keep' ) 
			? 'if'
			: ( $opts{mode} eq 'drop' ) 
				? 'unless'
				: do { &tmsg( "[Err]Unknown -mode [$opts{mode}]. Should be keep/drop."); die $!; }
		: 'if'
	; 
	$exe_code_src .= <<SRC; 
		$slct ( $DD_comp ) {
SRC
	$exe_code_src .= <<'SRC'; 
			print $out_fh $$fq_recR; 
			$num_keep ++; 
		} else {
			$num_drop ++; 
		}
		$tn ++; 
	}# End while 
SRC

	eval("$exe_code_src"); 

	&tmsg("[Rec] Total $tn reads, $num_keep kept + $num_drop dropped.\n"); 

	close ($src_fh); 
}#End for fname

&tmsg("[Rec] All over.\n"); 

sub readRecord {
	my $rec = readline($_[0]); 
	my $idx = ''; 
	$rec =~ m/\[(\d+)\]$/ and $idx = $1; 
	my $rdID = $rec; chomp($rdID); 
	$rdID =~ s/^\@// or die "[Err]Wrong read ID [$rdID]\n"; 
	$rec = $rec . readline($_[0]) . readline($_[0]) . readline($_[0]) ; 
	return(\$rec, $idx, $rdID); 
}

sub openFH ($$) {
	my $f = shift;
	my $type = shift;
	my %goodFileType = qw(
		<       read
		>       write
		read    read
		write   write
	);
	defined $type or $type = 'read';
	defined $goodFileType{$type} or die "[Err]Unknown open method tag [$type].\n";
	$type = $goodFileType{$type};
	local *FH;
	if ($type eq 'read') {
		if ($f =~ m/\.gz$/) {
			open (FH, '-|', "gzip -cd $f") or die "[Err]$! [$f]\n";
			# open ($tfh, '-|', "gzip -cd $f") or die "[Err]$! [$f]\n";
		} elsif ( $f =~ m/\.bz2$/ ) {
			open (FH, '-|', "bzip2 -cd $f") or die "[Err]$! [$f]\n";
		} else {
			open (FH, '<', "$f") or die "[Err]$! [$f]\n";                                                                                                                                                }
		} elsif ($type eq 'write') {
			if ($f =~ m/\.gz$/) {
				open (FH, '|-', "gzip - > $f") or die "[Err]$! [$f]\n";
			} elsif ( $f =~ m/\.bz2$/ ) {
				open (FH, '|-', "bzip2 - > $f") or die "[Err]$! [$f]\n";
			} else {
				open (FH, '>', "$f") or die "[Err]$! [$f]\n";
			}
		} else {
			# Something is wrong.
			die "[Err]Something is wrong here.\n";
		}
	return *FH;
}#End sub openFH

