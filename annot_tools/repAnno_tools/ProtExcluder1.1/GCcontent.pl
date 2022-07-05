#! /usr/bin/perl -w

$usage = "GCcotent.pl input-fasta-file \n";

# to caculate the number of each nucleotide in DNA sequences

if (@ARGV < 1) {die $usage;}

open(FA, "$ARGV[0]") || die $usage;

$seq = "";
$tctn = 0;
$tcta = 0;
$tctc = 0;
$tctg = 0;
$tctt = 0;

while (<FA>) {
    if (/>\s*(\S+)\s(.*)/) {
        $ctn = 0;
        $cta = 0;
        $ctc = 0;
        $ctg = 0;
        $ctt = 0;
	if ($seq) {
	    @sym = split(//, $seq);
	    foreach $sym (@sym) {
		if ($sym eq "T" || $sym eq "t") {
		    $ctt ++;
		}
		else {
		    if ($sym eq "G" || $sym eq "g") {
		$ctg ++;
	    }
		    else {
			if ($sym eq "C" || $sym eq "c") {
			    $ctc ++;
			}
			else {
			    if ($sym eq "A" || $sym eq "a") {
				$cta ++;
			    }
			    else {
				$ctn ++;
			    }
			}
		    }
		}
	    }
$total = $cta + $ctc + $ctg +$ctt;
	    if ($total > 0) {
$pergc = ($ctc +$ctg)/$total;	
}
	    else {
		$pergc = 0;
	    }
	printf "%06d %06d %06d %06d %06d %.3f\n",$cta,$ctc,$ctg,$ctt,$ctn,$pergc;
	}
	printf "%s\t", $1;
	$tctn = $tctn + $ctn;
	$tcta = $tcta + $cta;
	$tctc = $tctc + $ctc;
	$tctg = $tctg + $ctg;
	$tctt = $tctt + $ctt;
	$seq = "";
    }
     else {
	chomp;
	$seq .= $_;
    }
}
close FA;

@sym = split(//, $seq);
$ctn = 0;
$cta = 0;
$ctc = 0;
$ctg = 0;
$ctt = 0;
foreach $sym (@sym) {
	if ($sym eq "T" || $sym eq "t") {
	    $ctt ++;
	}
	else {
	    if ($sym eq "G" || $sym eq "g") {
                $ctg ++;
            }
	    else {
		if ($sym eq "C" || $sym eq "c") {
		    $ctc ++;
		}
		else {
		    if ($sym eq "A" || $sym eq "a") {
			$cta ++;
		    }
		    else {
			$ctn ++;
		    }
		}
	    }
	}
}

$total = $cta + $ctc + $ctg +$ctt;
if ($total > 0) {
$pergc = ($ctc +$ctg)/$total;	
printf "%06d %06d %06d %06d %06d %.3f\n",$cta,$ctc,$ctg,$ctt,$ctn,$pergc;
}
$tctn = $tctn + $ctn;
$tcta = $tcta + $cta;
$tctc = $tctc + $ctc;
$tctg = $tctg + $ctg;
$tctt = $tctt + $ctt;

$tc = $tcta +$tctc +$tctg +$tctt;
$tcn = $tctn + $tc;
print "A C G T N totalnoN total\n";

printf "%08d %08d %08d %08d %08d %08d %08d\n", $tcta,$tctc,$tctg,$tctt,$tctn,$tc,$tcn;

printf "AT %08d GC %08d\n", ($tcta+$tctt),($tctc+$tctg);
$gc = ( $tc == 0 ) ? 0 : ($tctg +$tctc)*100/$tc;
printf "GC is %.3f\n",$gc;
