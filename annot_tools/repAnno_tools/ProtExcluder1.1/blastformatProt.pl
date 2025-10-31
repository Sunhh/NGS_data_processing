#! /usr/bin/perl
# 20251031: HS edit to fit new BlastX output.

$usage = "blastformatProt.pl blastx_output_file score_cutoff identity_cutoff\nwhere the latter two are optional (default 0)\n";

# to reformat blastx ouput

if (@ARGV < 1) {die "$usage";}
if (@ARGV > 1) {$score_cutoff = $ARGV[1];}
else {$score_cutoff = 0;}
if (@ARGV > 2) {$iden_cutoff = $ARGV[2];}
else  {$iden_cutoff = 0;}
open(BLT, "$ARGV[0]") || die "Can not open BLAST output $ARGV[0].\n$usage";

$score = -1;
while (<BLT>) {
    if (/^Query=\s+(\S+)/) {
	$query_current = $1;
	$ct = "q";
    }
    elsif (/^>\s*(\S+)/) {
	$sbjct_current = $1;
	$ct = "sb";
    }
    elsif (/^Length=(\d+)/) {
	if ($ct eq "q") {
	    $q_len_current = $1;
	}
	elsif ($ct eq "sb") {
	    $sb_len_current = $1;
	}
    }
    elsif (/^\s+Score =\s*(\d+)/) {
	if (&filter) {
	    &report;
	}
	$head_q_boolean = 1;
	$head_s_boolean = 1;
	$query = $query_current;
	$q_len = $q_len_current;
	$sb_len = $sb_len_current;
	$sbjct = $sbjct_current;
	$score = $1;
    }
    elsif (/^\s+Identities =\s+(\d+)\/(\d+)/) {
	$iden = $1/$2*100;
    }
    elsif (/^Query\s+(\d+)\s+(\S+)\s+(\d+)/) {
	if ($head_q_boolean) {
	    $head_query = $1;
	    $head_q_boolean = 0;
	}
	$tail_query = $3;
    }
    elsif (/^Sbjct\s+(\d+)\s+(\S+)\s+(\d+)/) {
	if ($head_s_boolean) {
	    $head_sbjct = $1;
	    $head_s_boolean = 0;
	}
	$tail_sbjct = $3;
    }
}

close BLT;

if (&filter) {
    &report;
}

sub filter {
    if ($score > 0 && $score > $score_cutoff && $iden > $iden_cutoff) {
	return 1;
    }
    return 0;
}


sub report {
    if ($head_query < $tail_query) {
    printf("%06d %03d %05d %05d %-5d %-30s %05d %05d %05d %-30s +\n", $score, $iden, $head_query, 
$tail_query, $q_len, $query, $head_sbjct, $tail_sbjct, $sb_len, $sbjct);
    }
    else {
 printf("%06d %03d %05d %05d %-5d %-30s %05d %05d %05d %-30s C\n", $score, $iden, $tail_query,
	$head_query, $q_len, $query, $head_sbjct, $tail_sbjct, $sb_len, $sbjct);
    }

}
