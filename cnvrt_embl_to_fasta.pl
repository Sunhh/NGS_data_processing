#!/usr/bin/perl
use strict; 
-t and !@ARGV and die "perl $0 in.embl > out.fasta\n"; 

my $want_type='tab'; 
$want_type = 'fasta'; 
# >tr|AC|ID DE OS=OS GN=GN PE=PE SV=SV 
my @need_lab=qw/AC ID DE OS GN PE SV SQ/; 
my %info; 
for (@need_lab) {
	$info{$_} = ''; 
}

my $prev_lab = ''; 

while (<>) {
	chomp; 
	if (m!^//$!) {
		# Changing information; 
		if ( $info{'GN'} ne '' ) {
			$info{'GN'} =~ m!^(?:ORF|OrderedLocus|GLYI)?Name(?:s)?=([^;]+);! or die "GN:$info{GN}\n"; 
			$info{GN} = $1; 
		}
		if ($info{PE} ne '') {
			$info{PE} =~ m!^(\d+):! or die "PE:$info{PE}\n"; 
			$info{PE} = $1; 
		}
		if ($info{DE} ne '') {
			$info{DE} =~ s!^.*?Full=!! or die "DE:$info{DE}\n"; 
			my $flg = ''; 
			$info{DE} =~ m!Flags: ([^;]+);! and $flg = " ($1)"; 
			$info{DE} =~ s!\s*[^;\s]+\=.*$!!; 
			$info{DE} =~ s!;(?:SubName|Flags|AltName):.*$!!; 
			$info{DE} =~ s!;$!!; 
			$info{DE} .= $flg; 
		}
		if ($info{OS} ne '') {
			$info{OS} =~ s!\.$!!; 
			$info{OS} =~ s! \([^()]+\)!!g; 
		}
		if ($info{SQ} ne '') {
			$info{SQ} =~ s!\s!!g; 
		}
		# output and clear
		$prev_lab = ''; 
		if ($want_type eq 'tab') {
			for my $tl (@need_lab) {
				print STDOUT "$tl\t$info{$tl}\n"; 
				$info{$tl} = ''; 
			}
			print STDOUT "//\n"; 
		} elsif ($want_type eq 'fasta') {
			my $outTxt = ">tr|$info{AC}|$info{ID} $info{DE}"; 
			defined $info{OS} and $info{OS} ne '' and $outTxt .= " OS=$info{OS}"; 
			defined $info{GN} and $info{GN} ne '' and $outTxt .= " GN=$info{GN}"; 
			defined $info{PE} and $info{PE} ne '' and $outTxt .= " PE=$info{PE}"; 
			defined $info{SV} and $info{SV} ne '' and $outTxt .= " SV=$info{SV}"; 
			my $oS = $info{SQ}; 
			$oS =~ s/(.{60})/$1\n/g; chomp($oS); 
			print STDOUT "$outTxt\n$oS\n"; 
		}
		for my $tl (@need_lab) {
			$info{$tl} = ''; 
		}
	} elsif (s!^(\S{2})\s+!!) {
		my $tl = $1; 
		$prev_lab = $tl; 
		if ($tl eq 'ID') {
			m!^(\S+)! or die "ID:$_\n"; 
			$info{ID} = $1; 
		} elsif ($tl eq 'AC') {
			m!^(\S+);! or die "AC:$_\n"; 
			$info{AC} = $1; 
		} elsif ($tl eq 'DE') {
			$info{DE} .= $_; 
		} elsif ($tl eq 'OS') {
			$info{OS} ne '' and $info{OS} .= ' '; 
			$info{OS} .= $_; 
		} elsif ($tl eq 'GN') {
			$info{GN} .= $_; 
		} elsif ($tl eq 'PE') {
			$info{PE} .= $_; 
		} elsif ($tl eq 'DT') {
			if ($info{SV} eq '' and $_ =~ m!sequence version (\d+)\.!) {
				$info{SV} = $1; 
			} elsif ( $_ =~ m!sequence version (\d+)\.! ) {
				warn "SV repeat: $_\n"; 
			}
		} elsif ($tl eq 'SQ') {
			; 
		}
	} elsif (s!^\s+!!) {
		if ($prev_lab eq 'SQ') {
			$info{SQ} .= $_; 
		}
	} else {
		# Other types; 
		; 
	}
}
