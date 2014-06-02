#!/usr/bin/perl -w 
# edit 2010-09-24
# edit 2010-09-25 fit for blastn -query -subject command. 
# edit 2010-09-26 not to record sequence until -snp is defined. 
# edit 2010-10-07 using a 'index' method to do regexp match. 
# edit 2010-11-08 fix a bug occurring when there is a query without any sbjct hits. The line is like "^***** No hits found *****\n"; 
# edit 2014-06-02 Now can be used for both results from blast and blast+, and can give infor for subject_description. 

use strict; 

#-t and !@ARGV and die "perl $0 *.bn0\n"; 
-t and !@ARGV and &usage(); 

use Getopt::Long;
my %opts; 
GetOptions(\%opts,"help!","in:s","out:s","snp!","bn6:s"); 
sub usage {
	print STDOUT <<HELP; 
*************************** Instruction ***************************
Usage: perl $0 -in *.bn0 
-out    file name for bn6 output. 
-snp    calculate SNP positions if given. incompatible with -out. 
-bn6    bn6 file to restrict SNP finding. incompatible with -out. 
*************************** Instruction ***************************
HELP
	exit(0); 
}#end for usage(); 

sub openFH {
	my $f = shift;
	my $type = shift; (defined $type and $type =~ /^[<>|]+$/) or $type = '<';
	local *FH;
	open FH,$type,"$f" or die "Failed to open file [$f]:$!\n";
	return *FH;
}# end sub openFH

$opts{help} and &usage(); 

# region tag. 
my $region = 'head'; # qw/query sbjct head alignment/
my (%info, @ta, %snpaln); 
my $obn6fh = \*STDOUT; 

if ($opts{snp}) {
	open BN6, '<', "$opts{bn6}" or do { print "failed to open file(-bn6) $opts{bn6}. $!\n"; &usage(); }; 
	my %uniq; 
	while (<BN6>) {
		s/[^\S\t]+$//;
		/^\s*$/ and next; 
		my @tmp = split(/\t/, $_); 
		defined $uniq{$_} or push(@{$snpaln{$tmp[0]}}, $_); 
		$uniq{$_} = 1; 
	}
	close BN6; 
}else{
	$obn6fh = &openFH($opts{out}, '>'); 
	print $obn6fh join("\t", qw/qseqid sseqid pident aln_len mismatch gapopen qstart qend sstart send evalue bitscore qlen slen strand sdef/)."\n"; 
}


#    When not provided, the default value is:
#   'qseqid sseqid pident length mismatch gapopen qstart qend sstart send
#   evalue bitscore', which is equivalent to the keyword 'std'
open BN0,'<',"$opts{in}" or do { print "failed to open file(-in) $opts{in}. $!\n"; &usage(); }; 
while (<BN0>) {
	s/[^\S\t]+$//; 
	if ($region eq 'head') {
		if (@ta = /^Query\s*=\s*(\S+)\s*(.*)/i) {
			if (defined $info{send}) {
				defined $info{qstr} or $info{qstr} = 'NA'; defined $info{sstr} or $info{sstr} = 'NA'; 
				my $oline = join("\t", &outline(), @info{qw/qlen slen/}, "$info{qstr}$info{sstr}", $info{sdef}); 
				$opts{snp} or print $obn6fh "$oline\n"; 
				if ($opts{snp}) {
					for my $snpl (@{$snpaln{$info{qid}}}) {
						index($oline, $snpl) == -1 and next; 
						# $oline =~ /^$snpl/ or next; 
						&listSNP(\%info); 
					}
				}
				$info{send} = undef(); 
			}
			$info{qid} = $ta[0]; 
			$info{qdef} = $ta[1]; 
			$info{get_qlen} = 0; 
			$info{get_sdef} = 0; 
			$region = 'query'; 
		}
	}elsif ($region eq 'query') {
		#if (@ta = /^>\s*(\S+)\s*(.*)/) {
		if (@ta = /^>\s*(\S+)\s*(.*)/ or @ta = /^Subject\s*=\s*(\S+)\s*(.*)/i) {
			if (defined $info{send}) {
				defined $info{qstr} or $info{qstr} = 'NA'; defined $info{sstr} or $info{sstr} = 'NA'; 
				my $oline = join("\t", &outline(), @info{qw/qlen slen/}, "$info{qstr}$info{sstr}", $info{sdef}); 
				$opts{snp} or print $obn6fh "$oline\n"; 
				if ($opts{snp}) {
					for my $snpl (@{$snpaln{$info{qid}}}) {
						index($oline, $snpl) == -1 and next; 
						# $oline =~ /^$snpl/ or next; 
						&listSNP(\%info); 
					}
				}
				$info{send} = undef(); 
			}
			$info{sid} = $ta[0]; 
			$info{sdef} = $ta[1]; 
			$info{get_slen} = 0; 
			$info{get_sdef} = 0; 
			$region = 'sbjct'; 
		}elsif (@ta = /^Length=([\d,]+)/i or @ta= m/^\s{8,}\(([\d,]+) letters\)$/i) {
			$info{qlen} = $ta[0]; $info{qlen} =~ s/,//g; 
			$info{get_qlen} = 1; 
		}elsif ($info{get_qlen} == 0) {
			$info{qdef} .= " $_"; 
		}elsif (/^\*+/) { 
			# edit on 2010-11-08 
			# here the line is '***** No hits found *****'; 
			$region = 'head'; # return to head type. 
		}
	}elsif ($region eq 'sbjct') {
		if (@ta = /^Query\s*=\s*(\S+)\s*(.*)/i) {
			if (defined $info{send}) {
				defined $info{qstr} or $info{qstr} = 'NA'; defined $info{sstr} or $info{sstr} = 'NA'; 
				my $oline = join("\t", &outline(), @info{qw/qlen slen/}, "$info{qstr}$info{sstr}", $info{sdef}); 
				$opts{snp} or print $obn6fh "$oline\n"; 
				if ($opts{snp}) {
					for my $snpl (@{$snpaln{$info{qid}}}) {
						index($oline, $snpl) == -1 and next; 
						# $oline =~ /^$snpl/ or next; 
						&listSNP(\%info); 
					}
				}
				$info{send} = undef(); 
			}
			$info{qid} = $ta[0]; 
			$info{qdef} = $ta[1]; 
			$info{get_qlen} = 0; 
			$region = 'query'; 
		}elsif (@ta = />\s*(\S+)\s*(.*)/) {
			$info{sid} = $ta[0]; 
			$info{sdef} = $ta[1]; 
			$info{get_slen} = 0; 
			$info{get_sdef} = 0; 
		}elsif (@ta = /^\s*Length\s*=\s*([\d,]+)/i and $info{get_sdef} == 1) {
			$info{slen} = $ta[0]; $info{slen} =~ s/,//g; 
			$info{get_slen} = 1; 
		}elsif ($info{get_slen} == 0 and $info{get_sdef} == 0) {
			if (m/^$/) {
				$info{get_sdef} = 1; 
			}else{
				$_ =~ s/^\s+//; 
				$info{sdef} .= " $_"; 
			}
		}elsif (@ta = /^\s*Score += +(\S+) +bits +\((\d+)\), +Expect(?:\(\d+\))? += +(\S+)/i) {
			$info{score} = $ta[0]; 
			$info{bs} = $ta[1]; 
			($info{evalue} = $ta[2]) =~ s/,$//; 
			for my $ttt ( qw/qstart sstart qend send qseq sseq/ ) {
				$info{$ttt} = undef(); 
			}
		#}elsif (@ta = /^\s*Identities = (\d+)\/(\d+) \(([\d.]+)\%\), Gaps = (\d+)\/(\d+)/i) {
		}elsif ( @ta = m!^\s*Identities += +(\d+)\/(\d+) +\(([\d.]+)\%\)!i ) { 
#			@ta = /^\s*Identities += +(\d+)\/(\d+) +\(([\d.]+)\%\), +Gaps += +(\d+)\/(\d+)/i 
#			or @ta = /^\s*Identities += +(\d+)\/(\d+) +\(([\d.]+)\%\), +Positives += (\d+)\/(\d+) +\(([\d,]+)\%\), +Gaps += +(\d+)\/(\d+)/i
#			or @ta = m!^\s*Identities += +(\d+)\/(\d+) +\(([\d.]+)\%\)!i
#		) {  
			for my $ttt ( qw/match_base aln_len ident positives pos_ttl pos_ident gap_open gap_ttl/ ) {
				$info{$ttt} = undef(); 
			}
			@info{qw/match_base aln_len ident/} = @ta[0,1,2]; 
			if ( @ta = ($_ =~ m! Positives += (\d+)\/(\d+) +\(([\d,]+)\%\)!i) ) {
				@info{qw/positives pos_ttl pos_ident/} = @ta[0,1,2]; 
			}
			if ( @ta = ($_ =~ m! Gaps += +(\d+)\/(\d+)!i) ) {
				@info{qw/gap_open gap_ttl/} = @ta[0,1]; 
			}
			if (!defined $info{gap_open}) {
				$info{gap_open}=$info{gap_ttl}=0; 
			}
			if (!defined $info{positives}) {
				$info{pos_ttl} = $info{aln_len}; 
				$info{positives} = $info{match_base}; 
			}
			$info{mis_mat} = $info{aln_len}--$info{match_base}--$info{gap_open}; 
		}elsif (@ta = /^\s*Strand=(\w+)\/(\w+)/i or @ta = /^\s*Frame += +([\+\-\d]+)(?:\/([\+\-\d]+))?/i) {
			defined $ta[1] or $ta[1] = ''; 
			$info{qstr} = $ta[0]; $info{sstr} = $ta[1]; 
			$info{qstr} =~ s/Plus/+/i; $info{qstr} =~ s/Minus/\-/i; 
			$info{sstr} =~ s/Plus/+/i; $info{sstr} =~ s/Minus/\-/i; 
		}elsif (@ta = /^Query[\s:]+(\d+)\s+(\S+)\s+(\d+)/i) {
			defined $info{qstart} or $info{qstart} = $ta[0]; 
			$info{qseq} .= $ta[1]; 
			defined $info{qend} or $info{qend} = $ta[2]; 
			$info{qend} < $ta[2] and $info{qend} = $ta[2]; 
			$region = 'alignment'; 
		}else{
		}
	}elsif ($region eq 'alignment') {
		if (@ta = /^Query\s*=\s*(\S+)\s*(.*)/i) {
			if (defined $info{send}) {
				defined $info{qstr} or $info{qstr} = 'NA'; defined $info{sstr} or $info{sstr} = 'NA'; 
				my $oline = join("\t", &outline(), @info{qw/qlen slen/}, "$info{qstr}$info{sstr}", $info{sdef}); 
				$opts{snp} or print $obn6fh "$oline\n"; 
				if ($opts{snp}) {
					for my $snpl (@{$snpaln{$info{qid}}}) {
						index($oline, $snpl) == -1 and next; 
						# $oline =~ /^$snpl/ or next; 
						&listSNP(\%info); 
					}
				}
				$info{send} = undef(); 
			}
			$info{qid} = $ta[0]; 
			$info{qdef} = $ta[1]; 
			$info{get_qlen} = 0; 
			$region = 'query'; 
		#}elsif (@ta = /^>\s*(\S+)\s*(.*)/) {
		}elsif (@ta = /^>\s*(\S+)\s*(.*)/ or @ta = /^Subject\s*=\s*(\S+)\s*(.*)/i) {
			if (defined $info{send}) {
				defined $info{qstr} or $info{qstr} = 'NA'; defined $info{sstr} or $info{sstr} = 'NA'; 
				my $oline = join("\t", &outline(), @info{qw/qlen slen/}, "$info{qstr}$info{sstr}", $info{sdef}); 
				$opts{snp} or print $obn6fh "$oline\n"; 
				if ($opts{snp}) {
					for my $snpl (@{$snpaln{$info{qid}}}) {
						index($oline, $snpl) == -1 and next; 
						# $oline =~ /^$snpl/ or next; 
						&listSNP(\%info); 
					}
				}
				$info{send} = undef(); 
			}
			$info{sid} = $ta[0]; 
			$info{sdef} = $ta[1]; 
			$info{get_slen} = 0; 
			$info{get_sdef} = 0; 
			$region = 'sbjct'; 
		}elsif (@ta = /^\s*Score += +(\S+) +bits +\((\d+)\), +Expect(?:\(\d+\))? += +(\S+)/i) {
# editting here. 
			if (defined $info{send}) {
				defined $info{qstr} or $info{qstr} = 'NA'; defined $info{sstr} or $info{sstr} = 'NA'; 
				my $oline = join("\t", &outline(), @info{qw/qlen slen/}, "$info{qstr}$info{sstr}", $info{sdef}); 
				$opts{snp} or print $obn6fh "$oline\n"; 
				if ($opts{snp}) {
					for my $snpl (@{$snpaln{$info{qid}}}) {
						index($oline, $snpl) == -1 and next; 
						# $oline =~ m/^$snpl/ or next; 
						&listSNP(\%info); 
					}
				}
				$info{send} = undef(); 
			}
			$info{score} = $ta[0]; 
			$info{bs} = $ta[1]; 
			($info{evalue} = $ta[2]) =~ s/,$//; 
			$info{qstart} = undef(); $info{sstart} = undef(); 
			$info{qend} = undef(); $info{send} = undef(); 
			$info{qseq} = undef(); $info{sseq} = undef(); 
			$region = 'sbjct'; 
		}elsif (@ta = /^Query[\s:]+(\d+)\s+(\S+)\s+(\d+)/i) {
			defined $info{qstart} or $info{qstart} = $ta[0]; 
			$opts{snp} and $info{qseq} .= $ta[1]; 
			defined $info{qend} or $info{qend} = $ta[2]; 
			$info{qend} < $ta[2] and $info{qend} = $ta[2]; 
		}elsif (@ta = /^\s*Sbjct[\s:]+(\d+)\s+(\S+)\s+(\d+)/i) {
			defined $info{sstart} or $info{sstart} = $ta[0]; 
			$opts{snp} and $info{sseq} .= $ta[1]; 
			$info{send} = $ta[2]; 
			# defined $info{send} or $info{send} = $ta[2]; 
			# $info{send} > $ta[2] and $info{send} = $ta[2]; 
		}elsif (/^Lambda/) {
			if (defined $info{send}) {
				defined $info{qstr} or $info{qstr} = 'NA'; defined $info{sstr} or $info{sstr} = 'NA'; 
				my $oline = join("\t", &outline(), @info{qw/qlen slen/}, "$info{qstr}$info{sstr}", $info{sdef}); 
				$opts{snp} or print $obn6fh "$oline\n"; 
				if ($opts{snp}) {
					for my $snpl (@{$snpaln{$info{qid}}}) {
						index($oline, $snpl) == -1 and next; 
						# $oline =~ /^$snpl/ or next; 
						&listSNP(\%info); 
					}
				}
				%info = (); 
			}# for lambda 
		}else{
			; 
		}
	}
}#end while for BN0
close BN0; 
$opts{snp} or close ($obn6fh); 

sub outline {
	return join("\t", @info{qw/qid sid ident aln_len mis_mat gap_open qstart qend sstart send evalue score/}); 
}

# list all SNPs in the records; 
sub listSNP {
	my $rr = shift; 
	defined $rr->{qseq} or die "Wrong!\n$rr->{sid} have no qseq.\n"; 
	my @qs = split(//, $rr->{qseq}); 
	my @ss = split(//, $rr->{sseq}); 
	my ($q_delt, $s_delt) = (0, 0); 
	for (my $i=0; $i<@qs; $i++) {
		defined $ss[$i] or die "Failed to find [$i+1] base for sseq.\n"; 
		$qs[$i] eq '-' or $q_delt++; 
		$ss[$i] eq '-' or $s_delt++; 
		if ($qs[$i] !~ /$ss[$i]/i) {
			my ($qpos, $spos); 
			$qpos = $rr->{qstart} + $q_delt-1; 
			if ($rr->{sstr} eq '+' or ($rr->{sstr} eq 'NA' and $rr->{qstr} eq 'NA')) {
				$spos = $rr->{sstart}+$s_delt-1; 
			}else{
				$spos = $rr->{sstart}-$s_delt+1; 
			}
			print STDOUT join("\t", $qpos, $qs[$i], $ss[$i], $spos, @{$rr}{qw/qid sid qlen slen qstr sstr sdef/})."\n"; 
		}
	}
}# end sub listSNP; 




