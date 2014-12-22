#!/usr/bin/perl
use strict; 
use warnings; 
use File::Basename; 
use LogInforSunhh; 

-t and !@ARGV and &stopErr("[Err] perl $0 assembly_chunk0.err\n"); 

my %infor; 
my @ids = qw/asmDir bad_proc start_time end_time fill_len fill_shrt no_extender no_tails/; 
{for my $tk (@ids) {
	$infor{$tk} = $tk; 
}
&outputInfor(\%infor); 
%infor = (); 
}

my $is_problem = 0; 
while (<>) {
	chomp; 
	if (s/^(\d\d\d\d\-\d\d-\d+ \d+:\d+:\d+)\,(\d+)\s*//) {
		# This is normal records. 
		$is_problem = 0; 
		my ($t_time, $t_num) = ($1, $2); 
		if (m/^\[INFO\] Running \S+Assembly\.py (\S+)/) {
			my $t_asmDir = $1; 
			$t_asmDir = basename($t_asmDir); 
			if (keys %infor > 2) {
				&outputInfor(\%infor); 
			} elsif (keys %infor == 2) {
				if ( $infor{asmDir} ne $t_asmDir ) {
					&outputInfor(\%infor); 
					die "$infor{asmDir} ne $t_asmDir\n$_\n"; 
				} else {
					next; 
				}
			} elsif (keys %infor == 0) {
				# 
			} else {
				&stopErr("[Err] Why here! $_\n"); 
			}
			%infor = (); 
			$infor{asmDir} = $t_asmDir; 
			$infor{start_time} = $t_time; 
		} elsif (m/^\[INFO\] Finished/) {
			$infor{end_time} = $t_time; 
			&outputInfor(\%infor); 
			%infor = (); 
		} elsif (m/^\[INFO\] estimated fill len (\d+)/) { 
			$infor{fill_len} = $1; 
		} elsif (m/^\[WARNING\] (fill|f1e|f2e) sequence is small \((\d+)bp\) can\'t call consensus/) { 
			$infor{fill_shrt} = "$1:$2"; 
		} elsif (m/^\[INFO\] no overlap between extenders/) { 
			$infor{no_extender} = 1; 
		} elsif (m/^\[INFO\] No tails/) {
			$infor{no_tails} = 1; 
		} else {
			if (m/^\[INFO\]\s+(?:Extracting tails|Preunited |Parsed |Found |\d+ reads had double tails|Mapping Tails|\[0\, \'\[INFO\] |\d+ tails mapped|\d+ reads span|\d+ reads extend flank|estimated flank|Cleaning \d+ temp files|Consolidating alignments|fill seq is too short to call consensus|no mapping\.+\s+picking \S+ seq)/) {
				; 
			} else {
				&stopErr("[Err] Check this line: $_\n"); 
			}
		}
	} elsif (m/^Traceback \(most recent call last\)/) {
		$infor{bad_proc} = 1; 
		$is_problem = 1; 
	} elsif ($is_problem == 1) { 
		next; 
	} else {
		&stopErr("[Err] Why here? : $_\n"); 
	}
}

# Input : \%infor : {@ids}
# asmDir     : ref0105227e5_ref0293445e3
# start_time : 2014-12-17 10:01:17
# end_time   : 2014-12-17 10:01:21
# fill_len   : 2005
# fill_shrt  : 0/1
# bad_proc   : 0/1
# no_extender : 0/1
# no_tails   : 0/1 
sub outputInfor {
	my $t_hash = shift; 
	my $back_line; 
	for my $tk (@ids) {
		defined $t_hash->{$tk} or $t_hash->{$tk} = -1; 
	}
	print STDOUT join("\t", @{$t_hash}{@ids})."\n"; 
}#End sub outputInfor 

