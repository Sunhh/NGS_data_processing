#!/usr/bin/perl
use strict; 
use warnings; 
use Statistics::Descriptive::Weighted; 

-t and !@ARGV and &stopErr("samtools view -h in.bam | perl $0 \n"); 


my $min_scf_len = 20e3; 
$min_scf_len = 1; 
$min_scf_len = 10e3; 

my $max_ins = 50e3; 

###################################################################
# BitPos  Flag  Chr Description
# 0 0x0001  p the read is paired in sequencing
# 1 0x0002  P the read is mapped in a proper pair
# 2 0x0004  u the query sequence itself is unmapped
# 3 0x0008  U the mate is unmapped
# 4 0x0010  r strand of the query (1 for reverse)
# 5 0x0020  R strand of the mate
# 6 0x0040  1 the read is the first read in a pair
# 7 0x0080  2 the read is the second read in a pair
# 8 0x0100  s the alignment is not primary
# 9 0x0200  f the read fails platform/vendor quality checks
# 10  0x0400  d the read is either a PCR or an optical duplicate
###### FLAG list prepared.

###### Make FLAG list
my %flag;
for (0 .. 2047) {
	my $binum = unpack("B32", pack("N", $_));
	$flag{$_} = [ reverse( split(//, sprintf("%011d", $binum) ) ) ];
}

my $opts_keep = '0=1,2=0,3=0,4=0,5=1'; 

my %use_flag; 

my %keep_flag; 
for ( sort { $a<=>$b } keys %flag ) {
	my @ta = @{$flag{$_}}; 
	CHK_TTL_KEEP:
	for my $exp1 ( split(/;/, $opts_keep) ) {
		my $should_keep = 1;
		CHK_EACH_KEEP:
		for my $exp2 ( split(/,/, $exp1) ) {
			$exp2 =~ m/^(\d+)=([01])$/ or die "Failed to parse |$exp2|\n";
			my ($coln, $colv) = ($1, $2);
			$ta[$coln] != $colv and do { $should_keep = 0; last CHK_EACH_KEEP; };
		}#End for my $exp2 [CHK_EACH_KEEP]

		# logical OR, so any of TRUE will make the keeping [CHK_TTL_KEEP] TRUE.
		$should_keep == 1 and do { $keep_flag{$_} = 1; last CHK_TTL_KEEP; };

	}#End for my $exp1 [CHK_TTL_KEEP]

	if ( defined $keep_flag{$_} and $keep_flag{$_} == 1 ) {
		defined $use_flag{$_} and &stopErr("[Err]Repeat flag[$_]\n"); 
		$use_flag{$_} = 1; 
	} else {
		$use_flag{$_} = 0; 
	}
}

print STDOUT join("\t", qw/ScaffID ScaffPos FragCov FragLenMean/)."\n"; 

my @frag_covs; 
my %scf_len; 
my $prev_scfID; 
my %scf_range; 
while (<>) {
	$. % 1e5 == 1 and &tsmsg("[Msg]Processing $. line.\n"); 
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($ta[0] =~ m/^\@/) {
		if ($ta[0] eq '@SQ') {
			$ta[1] =~ m/^SN:(\S+)$/ or &stopErr("[Err]Failed to parse1 $_\n"); 
			my $scfID = $1; 
			$ta[2] =~ m/^LN:(\d+)$/ or &stopErr("[Err]Failed to parse2 $_\n"); 
			my $scfLen = $1; 
			defined $scf_len{$scfID} and &stopErr("[Err]Repeat scafID $_\n"); 
			$scfLen >= $min_scf_len or next; 
			$scf_len{$scfID} = $scfLen ; 
		}
		next; 
	}
	my ($cur_flag, $scfID, $scfPosS, $mate_scf, $spanSize) = @ta[1,2,3,6,8]; 

	$use_flag{$cur_flag} == 1 or next; 
	$mate_scf eq '=' or next; 
	abs($spanSize) <= 50e3 or next; 
#	$spanSize > 0 or next; 
	my $scfPosE = $scfPosS + $spanSize - 1; 
	$scfPosE < $scfPosS and ($scfPosS, $scfPosE) = ($scfPosE, $scfPosS); 
	defined $scf_len{$scfID} or do { &tsmsg("[Msg] Skip short scaffold [$scfID]\n"); next; }; 
	if (defined $prev_scfID) {
		if ($prev_scfID eq $scfID) {
			for (my $p=$scfPosS; $p<=$scfPosE; $p++) {
				push(@{ $frag_covs[$p-1] }, $spanSize); 
				defined $scf_range{max} or $scf_range{max} = $p; 
				defined $scf_range{min} or $scf_range{min} = $p; 
				$scf_range{min} > $p and $scf_range{min} = $p; 
				$scf_range{max} < $p and $scf_range{max} = $p; 
			}
		} else {
			&tsmsg("[Msg]Output $prev_scfID len=$scf_len{$prev_scfID}\n"); 
			# for (my $p=0; $p<$scf_len{$prev_scfID}; $p++) {
			for (my $p=$scf_range{min}-1; $p<$scf_range{max}; $p++) {
				if ( scalar( @{$frag_covs[$p]} ) > 0 ) {
					print STDOUT join("\t", $prev_scfID, $p+1, 0, -1, '')."\n"; 
				} else {
					print STDOUT join( "\t", $prev_scfID, $p+1, scalar( @{$frag_covs[$p]} ), &mean( @{$frag_covs[$p]} ), join(';', @{$frag_covs[$p]}) )."\n"; 
				}
			}
			
			$prev_scfID = $scfID; 
			# $#frag_covs = $scf_len{$scfID}-1; 
			&tsmsg("[Msg]Initializing $scfID\n"); 
			for (my $p=0; $p<$scf_len{$scfID}; $p++) {
				$frag_covs[$p] = []; 
			}
			for (my $p=$scfPosS; $p<=$scfPosE; $p++) {
				push( @{$frag_covs[$p-1]}, $spanSize ); 
				defined $scf_range{max} or $scf_range{max} = $p; 
				defined $scf_range{min} or $scf_range{min} = $p; 
				$scf_range{min} > $p and $scf_range{min} = $p; 
				$scf_range{max} < $p and $scf_range{max} = $p; 
			}
		}
	}else{
		$prev_scfID = $scfID; 
		&tsmsg("[Msg]Initializing $scfID\n"); 
		for (my $p=0; $p<$scf_len{$scfID}; $p++) {
			$frag_covs[$p] = []; 
		}
		for (my $p=$scfPosS; $p<=$scfPosE; $p++) {
			push(@{$frag_covs[$p-1]}, $spanSize); 
			defined $scf_range{max} or $scf_range{max} = $p; 
			defined $scf_range{min} or $scf_range{min} = $p; 
			$scf_range{min} > $p and $scf_range{min} = $p; 
			$scf_range{max} < $p and $scf_range{max} = $p; 
		}
	}
}

# for (my $p=0; $p<$scf_len{$prev_scfID}; $p++) {
for (my $p=$scf_range{min}-1; $p<$scf_range{max}; $p++) {
	print STDOUT join("\t", $prev_scfID, $p+1, scalar( @{$frag_covs[$p]} ), &mean( @{$frag_covs[$p]} ), join(';', @{$frag_covs[$p]}) )."\n"; 
}

&tsmsg("[Rec]All done.\n"); 

sub mean {
	my $ttl = 0; 
	my $nn = 0; 
	scalar(@_) == 0 and return -1; 
	for (@_) {
		$_ >= $max_ins and next; 
		$ttl += $_; 
		$nn ++; 
	}
	return $ttl/$nn ; 
}

sub tsmsg {
	my $tt = scalar(localtime()); 
	print STDERR join('', "[$tt]", @_); 
}

sub stopErr {
	&tsmsg(@_); 
	exit(1); 
}

