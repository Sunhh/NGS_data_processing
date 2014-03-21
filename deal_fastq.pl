#!/usr/bin/perl 
# 2014-03-20 A script to deal with fastq format reads. It will be always in processing. 
#
use strict; 
use warnings; 
use List::Util qw(first max maxstr min minstr reduce shuffle sum); 
use Getopt::Long; 
my %opts; 

sub usage {
	print STDOUT <<HELP; 
#******* Instruction of this program *********#

Introduction:Deal with fasta format file. 
perl $0 in.fastq 

-fq2fa        [Boolean] Transformat fastq sequences to fasta format. 
-oQfile       [String] Name of quality file to output for -fq2fa fastq sequence. 

-fq2Val       [Boolean] Change the quality string to a set of quality numbers seprated by white space. 

-keep_len     [String] "min_len-max_len". Extract sequences whose lengths are between min_len and max_len.
-paired       [Boolean] Treat files as paired, and output reads in order of R1/R2 ... 

-sepR12_pref  [String] Input reads as R1-R2 order, and we output R1 to pref_R1.fq and R2 to pref_R2.fq
-joinR12      [String] Input files as paired, and join R1/R2 reads in a same stream. You can use \"-keep_len 0- -paired\" to get the same result. 

-showQscale   [Boolean] Show fastq quality scale value [33/64]. 

-33to64       [Boolean] Transform phred33 to phred64 quality. 
-64to33       [Boolean] Transform phred64 to phred33 quality. 

-rd_Num       [Boolean] Summary read number/length of fq files. 

-frag         [String] [start-end] output a fragment of the squence in single sequence fasta file;Start position is 1; May be s1-e1:s2-e2...
-frag_r       [Boolean] give out reversed string; 
-frag_c       [Boolean] give out complemented string as DNA; 

-phred_scale [INT] 33/64. Set scale value to get phred quality score. 

#******* Instruction of this program *********#
HELP
	exit(1); 
}

GetOptions(\%opts, 
	"phred_scale:i", 
	"fq2fa!", "oQfile:s", 
	"fq2Val!", 
	"keep_len:s", "paired!", 
	"sepR12_pref:s", 
	"joinR12!", # In fact, this function can be performed by keep_len() subroutine with "-keep_len 0- -paired "
	"showQscale!", 
	"33to64!", "64to33!", 
	"rd_Num!", 
	"frag:s", "frag_r!", "frag_c!", 
	"help!", 
); 

#****************************************************************#
#--------------Main-----Function-----Start-----------------------#
#****************************************************************#

&usage() if ( $opts{help} or ( -t and !@ARGV ) );


# Making File handles for reading;
our @InFp = () ; # 2007-8-29 16:07 È«¾Ö±äÁ¿!
if ( !@ARGV )
{
	@InFp = (\*STDIN);
}
else
{
	for (@ARGV) {
		push( @InFp, &openFH($_,'<') );
	}
}


## Global settings


## Call functions. 
&fq2fa() if ( $opts{fq2fa} ); 
&keep_len( $opts{keep_len} ) if ( defined $opts{keep_len} and $opts{keep_len} ne '' ); 
&showQscale() if ( $opts{showQscale} ); 
&fq2fq( 64-33 ) if ( $opts{'33to64'} ); 
&fq2fq( 33-64 ) if ( $opts{'64to33'} ); 
&sepR12( $opts{sepR12_pref} ) if ( defined $opts{sepR12_pref} and $opts{sepR12_pref} ne '' ); 
&joinR12() if ( $opts{joinR12} ); 
&fq2Val() if ( $opts{fq2Val} ); 
&rd_Num() if ( $opts{rd_Num} ); 
&fragmentRd() if ( defined $opts{frag} and $opts{frag} ne '' ); 

#****************************************************************#
#--------------Subprogram------------Start-----------------------#
#****************************************************************#

sub fragmentRd{
	my $fragR = &frag_regions($opts{frag}); 
	for my $fh ( @InFp ) {
		while ( my $rdRec = &get_fq_record($fh) ) {
			$rdRec->{seq} =~ s/\s//g; 
			$rdRec->{qual} =~ s/\s//g; 
			my $l_seq = length($rdRec->{seq}); 
			my ($strSeq, $strQual) = ('', ''); 
			my @Range; 
			for my $seR ( @$fragR ) {
				my ($add_s, $add_e) = @$seR; 
				$add_s < 0 and $add_s = $l_seq+$add_s+1; 
				$add_e eq 'end' and $add_e = $l_seq; 
				$add_e < 0 and $add_e = $l_seq+$add_e+1; 
				$add_e > $l_seq and $add_e = $l_seq; 
				if ( $add_s > $l_seq ) {
					; # Nothing happening. 
				} else {
					$strSeq  .= substr($rdRec->{seq},  $add_s-1, $add_e-$add_s+1); 
					$strQual .= substr($rdRec->{qual}, $add_s-1, $add_e-$add_s+1); 
				}
				push(@Range, "$add_s\-$add_e"); 
			}
			my $range = join(',', @Range); 
			if ( $opts{frag_c} ) { &rcSeq(\$strSeq, 'c'); $range = "C$range"; } 
			if ( $opts{frag_r} ) { $strSeq = reverse($strSeq); $strQual = reverse($strQual); $range = "R$range"; } 
			chomp($rdRec->{id}); $rdRec->{id} .= " [$range]"; 
			print STDOUT "\@$rdRec->{id}\n$strSeq\n+\n$strQual\n"; 
		}#while
	}#for
}#sub fragmentRd() 

sub showQscale {
	print STDOUT join("\t", qw/InFqName QualScale/)."\n"; 
	if ( !@ARGV ) {
		my $qscale = &guessPhredScale($InFp[0]); 
		print STDOUT join("\t", '[STDIN]', $qscale)."\n"; 
	} else {
		for (my $i=0; $i<@ARGV; $i++) {
			my $qscale = &guessPhredScale($InFp[$i]); 
			print STDOUT join("\t", $ARGV[$i], $qscale)."\n"; 
		}
	}
}#sub showQscale() 

sub fq2fa {
	my $oQfh = undef(); 
	if (defined $opts{oQfile} and $opts{oQfile} ne '') {
		$oQfh = &openFH($opts{oQfile}, '>'); 
	}
	my $rdNum = 0; 
	for my $fh ( @InFp ) {
		my $qscale = &guessPhredScale($fh); 
		if ( defined $oQfh ) {
			while ( my $rdRec = &get_fq_record($fh) ) {
				$rdNum ++; 
				$rdNum % 10e6 == 1 and &tsmsg("[Msg] Dealing $rdNum reads.\n"); 
				print STDOUT ">$rdRec->{id}$rdRec->{seq}"; 
				my $qVR = &qChar2qValue(\$rdRec->{qual}, $qscale); 
				print {$oQfh} ">$rdRec->{id}@$qVR\n"; 
			}
		} else {
			while ( my $rdRec = &get_fq_record($fh) ) {
				$rdNum ++; 
				$rdNum % 10e6 == 1 and &tsmsg("[Msg] Dealing $rdNum reads.\n"); 
				print STDOUT ">$rdRec->{id}$rdRec->{seq}"; 
			}
		}
	}# End for my $fh
	&tsmsg("[Rec] Finish $rdNum reads.\n"); 
}#sub fq2fa

sub rd_Num {
	print STDOUT join("\t", qw/InFile Total_size Total_Rd_num Mean_Rd_size Range_Rd_size PhredCut Time/)."\n"; 
	for (my $i=0; $i<@InFp; $i++) {
		my $fh = $InFp[$i]; 
		my $fn = ( defined $ARGV[$i] ) ? $ARGV[$i] : '[STDIN]' ; 
		my $qscale = &guessPhredScale($fh); 
		&tsmsg("[Rec] Reading file [$fn] with phred scale [$qscale]\n"); 
		my %dd; 
		$dd{rdN} = $dd{bpN} = 0; 
		$dd{len_range} = ['NA', 'NA']; 
		while ( my $rdRec = &get_fq_record($fh) ) {
			$dd{rdN} ++; $dd{rdN} % 10e6 == 1 and &tsmsg("[Msg] $dd{rdN} reads.\n"); 
			$rdRec->{seq} =~ s/\s//g; 
			my $tt_len = length( $rdRec->{seq} ); 
			$dd{bpN} += $tt_len; 
			if ($dd{len_range}[0] eq 'NA') {
				$dd{len_range}[0] = $tt_len; 
				$dd{len_range}[1] = $tt_len; 
			} else {
				if      ( $dd{len_range}[0] > $tt_len ) {
					$dd{len_range}[0] = $tt_len; 
				} elsif ( $dd{len_range}[1] < $tt_len ) {
					$dd{len_range}[1] = $tt_len; 
				}
			}
		}#while
		$dd{avgL} = ( $dd{rdN} == 0 ) ? 'NA' : $dd{bpN}/$dd{rdN} ; 
		print STDOUT join("\t", $fn, $dd{bpN}, $dd{rdN}, $dd{avgL}, "$dd{len_range}[0]-$dd{len_range}[1]", "Phred$qscale", scalar(localtime()))."\n"; 
		&tsmsg("[Rec] Finish $dd{rdN} reads in file [$fn].\n"); 
	}#for 
}#sub rd_Num()



sub fq2Val {
	local $" = " "; 
	for my $fh ( @InFp ) {
		my $qscale = &guessPhredScale($fh); 
		while ( my $rdRec = &get_fq_record($fh) ) {
			my $qVR = &qChar2qValue(\$rdRec->{qual}, $qscale); 
			print STDOUT "\@$rdRec->{id}$rdRec->{seq}+\n@$qVR\n"; 
		}
	}
}#sub fq2Val() 

# Filter reads by read length. 
sub keep_len {
	my ($min_len, $max_len) = &max_min_range( shift ); 

	if ( $opts{paired} ) {
		for (my $i=0; $i<@InFp; $i+=2) {
			my $fh1 = $InFp[$i]; 
			my $fh2 = $InFp[$i+1]; 
			RD: 
			while ( !eof($fh1) and !eof($fh2) ) {
				my $rdRec1 = &get_fq_record($fh1); 
				my $rdRec2 = &get_fq_record($fh2); 
				(my $ss1 = $rdRec1->{seq}) =~ s/\s//g; 
				my $ll1 = length($ss1); 
				$min_len > 0 and $ll1 < $min_len and next RD; 
				$max_len > 0 and $ll1 > $max_len and next RD; 
				(my $ss2 = $rdRec2->{seq}) =~ s/\s//g; 
				my $ll2 = length($ss2); 
				$min_len > 0 and $ll2 < $min_len and next RD; 
				$max_len > 0 and $ll2 > $max_len and next RD; 
				print STDOUT "\@$rdRec1->{id}$rdRec1->{seq}+\n$rdRec1->{qual}\n"; 
				print STDOUT "\@$rdRec2->{id}$rdRec2->{seq}+\n$rdRec2->{qual}\n"; 
			}#End while() RD: 
		}#End for (my $i=0; ...)
	} else {
		for my $fh (@InFp) {
			while ( my $rdRec = &get_fq_record($fh) ) {
				(my $ss = $rdRec->{seq}) =~ s/\s//g; 
				my $ll = length($ss); 
				$min_len > 0 and $ll < $min_len and next; 
				$max_len > 0 and $ll > $max_len and next; 
				print STDOUT "\@$rdRec->{id}$rdRec->{seq}+\n$rdRec->{qual}\n"; 
			}
		}#End for my $fh
	}# if ( $opts{paired} )
}# sub keep_len

sub joinR12 {
	for (my $i=0; $i<@InFp; $i+=2) {
		my $fh1 = $InFp[$i]; 
		my $fh2 = $InFp[$i+1]; 
		RD: 
		while ( !eof($fh1) and !eof($fh2) ) {
			my $rdRec1 = &get_fq_record($fh1); 
			my $rdRec2 = &get_fq_record($fh2); 
			print STDOUT "\@$rdRec1->{id}$rdRec1->{seq}+\n$rdRec1->{qual}\n"; 
			print STDOUT "\@$rdRec2->{id}$rdRec2->{seq}+\n$rdRec2->{qual}\n"; 
		}#End while() RD: 
	}#End for
}#sub joinR12() 



sub sepR12 {
	my $pref = shift; 
	my $of1h = &openFH( "${pref}_R1.fq", '>' ); 
	my $of2h = &openFH( "${pref}_R2.fq", '>' ); 
	for my $fh ( @InFp ) {
		while ( !eof($fh) ) {
			my $rdRec1 = &get_fq_record($fh); 
			my $rdRec2 = &get_fq_record($fh); 
			print {$of1h} "\@$rdRec1->{id}$rdRec1->{seq}+\n$rdRec1->{qual}\n"; 
			print {$of2h} "\@$rdRec2->{id}$rdRec2->{seq}+\n$rdRec2->{qual}\n"; 
		}# while ( !eof() )
	}# for my $fh 
}#sub sepR12() 


sub fq2fq {
	my $scale = shift; 
	defined $scale or $scale = 33-64; 
	for my $fh ( @InFp ) {
		while ( my $rdRec = &get_fq_record($fh) ) {
			my $new_qStrR = &char2char( \$rdRec->{qual}, $scale ); 
			print STDOUT "\@$rdRec->{id}$rdRec->{seq}+\n$$new_qStrR\n"; 
		}
	}#
}#sub fq2fq() 


#****************************************************************#
#--------------InnerSubprogram------------Start------------------#
#****************************************************************#

# phred33 to phred64 / phred64 to phred33
# Input : Reference of quality string, and the scale to change char to char. Default is to change phred64 to phred33. 
# Return: Reference of new quality string. 
sub char2char {
	my $raw_qR = shift; 
	my $scale = shift; 
	defined $scale or $scale = 33-64; 
	my @new_ascii; 
	for my $tv ( unpack("C*", $$raw_qR) ) {
		push( @new_ascii, $tv+$scale ); 
	}
	my $back_q = pack("C*", @new_ascii); 
	return(\$back_q); 
}#sub char2char() 


# transform illumina quality character string to phred values. 
# Input : a reference of quality string. And a phred scale used to cut ord() number. 
# Return: an array reference [q1,q2,q3,q4,...]
sub qChar2qValue {
	my $cR = shift; 
	my $scale = shift; 
	defined $scale or $scale = 33; 
	my @back = unpack("C*", $$cR); 
	for my $tv ( @back ) {
		$tv -= $scale; 
	}# for my $tv (@back)
#	pos($$cR) = 0; 
#	while ( $$cR =~ m/(\S)/g ) {
#		push(@back, ord($1)-$scale); 
#	}
	return(\@back); 
}#sub qChar2qValue() 

# Input : an array reference of quality values [q1,q2,q3,...]. And a phred scale used to get character. 
# Return: a reference of character string for the input values. 
sub qValue2qChar {
	my $vR = shift; 
	my $scale = shift; 
	defined $scale or $scale = 33; 
	my $back; 
	for my $tv (@$vR) {
		$back .= chr( $tv + $scale ); 
	}
	return (\$back); 
}

# Read one fastq from file handle once a time. 
# Input : a fastq file handle. 
# Return: A hash reference (\%) with keys: qw/id seq qual/, in which the first '@' in "id" is trimmed, and the quality string is chomped. 
sub get_fq_record {
	my $fh = shift; 
	my %back; 
	eof($fh) and return( undef() ); 
	$back{id} = readline($fh); 
	$back{id} =~ s/^\@// or &stopErr("[Err] The read id line seems wrong:\t$back{id}\n"); 
	$back{seq} = readline($fh); 
	readline($fh); 
	$back{qual} = readline($fh); 
	chomp($back{qual});
	return(\%back); 
}

sub tsmsg {
	my $tt = scalar(localtime()); 
	print STDERR join('', "[$tt] ", @_); 
}

sub stopErr {
	&tsmsg(@_); 
	exit(1); 
}

# Guess phred scale from file handle input. 
# Input : File handle of fastq file. 
# Return: Phred scale of fastq quality format. 
sub guessPhredScale {
	defined $opts{phred_scale} and return ($opts{phred_scale}); 
	my $fh = shift; 
	my $not_get = 1; 
	while ($not_get) {
		my $rd = &get_fq_record($fh); 
		chomp($rd->{qual}); 
		pos($rd->{qual}) = 0; 
		while ( $rd->{qual} =~ m/(\S)/g ) {
			my $tQV = ord($1); 
			if ( $tQV <= 58 ) {
				seek( $fh, 0, 0 ); 
				return ( 33 ); # -phred33
			} elsif ( $tQV >= 75 ) {
				seek( $fh, 0, 0 ); 
				return ( 64 ); # -phred64
			}
		}
	}
	&tsmsg("[Err] Failed to guess phred scale for file handle [$fh], use 33 as default.\n"); 
	seek( $fh, 0, 0 ); 
	return( 33 ); # -phred33
}


# input ($seq_ref, $deal_tag); deal_tag : 'r' => reverse, 'c' => complemented, 'rc' => reverse and complemented; Default 'rc';
# no output, edit the input sequence reference.
sub rcSeq {
	my $seq_r = shift;
	my $tag = shift; defined $tag or $tag = 'rc'; # $tag = lc($tag);
	my ($Is_r, $Is_c) = (0)x2;
	$tag =~ /r/i and $Is_r = 1;
	$tag =~ /c/i and $Is_c = 1;
	!$Is_r and !$Is_c and die "Wrong Input for function rcSeq! $!\n";
	$Is_r and $$seq_r = reverse ($$seq_r);
	$Is_c and $$seq_r =~ tr/acgturykmbvdhACGTURYKMBVDHwWsSnN/tgcaayrmkvbhdTGCAAYRMKVBHDwWsSnN/; # edit on 2013-09-11 No difference in result.
	return 0;
}# 2007-9-11 9:46 ÖÆ×÷¶ÔÓ¦·´Ïò»¥²¹ÐòÁÐ;
#        a       a; adenine
#        c       c; cytosine
#        g       g; guanine
#        t       t; thymine in DNA; uracil in RNA
#        m       a or c
#        k       g or t
#        r       a or g
#        y       c or t
#        w       a or t
#        s       c or g
#        v       a or c or g; not t
#        b       c or g or t; not a
#        h       a or c or t; not g
#        d       a or g or t; not c
#        n       a or c or g or t

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
	# my $tfh;
	if ($type eq 'read') {
		if ($f =~ m/\.gz$/) {
			open (FH, '-|', "gzip -cd $f") or die "[Err]$! [$f]\n";
			# open ($tfh, '-|', "gzip -cd $f") or die "[Err]$! [$f]\n";
		} elsif ( $f =~ m/\.bz2$/ ) {
			open (FH, '-|', "bzip2 -cd $f") or die "[Err]$! [$f]\n";
		} else {
			open (FH, '<', "$f") or die "[Err]$! [$f]\n";
		}
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

# Input : "min_len-max_len"
# Return: (min_len, max_len)
sub max_min_range {
	my $kplen = shift; 
	my ($min_len, $max_len) = (-1, -1); 
	if ( $kplen =~ m!^\s*(\d+)\-(\d+)\s*$! ) {
		($min_len, $max_len) = ($1, $2);
	} elsif ( $kplen =~ m!^\s*\-(\d+)\s*$! ) {
		$max_len = $1;
	} elsif ( $kplen =~ m!^\s*(\d+)\-\s*$! ) {
		$min_len = $1;
	} else {
		&stopErr( "[Err] Failed to parse option [-keep_len $kplen]\n" );
	}
	return($min_len, $max_len); 
}# sub max_min_range

# Input : string to give a set of ranges. 
# Return: A reference of array whose element is array of [start, end]
sub frag_regions {
	my $frag_str = shift; 
	( defined $frag_str and $frag_str ne '') or &stopErr("[Err] No value given to \$frag_str.\n"); 
	my @back; 
	for ( split(/:/, $frag_str) ) {
		s/^\s+//; s/\s+$//; 
		if ( m/^(\-?(?:\d+)?)-(\-?(?:\d+)?)$/ ) {
			my ($ss, $ee) = ($1,$2); 
			push( @back, [ ($ss) ? $ss : 1, ($ee) ? $ee : 'end' ] ); 
			(defined $ss and defined $ee and $ss > $ee and $ee > 0) and &stopErr("[Err] Bad fragment region [$frag_str].\n"); 
		} else {
			&stopErr( "[Err] fragment region [$frag_str] not accepted.\n" ); 
		}
	}
	return(\@back); 
}#sub frag_regions() 


