#! /usr/bin/perl -w

$usage = "countaanu.pl blastx.mt file\n";

# to calculate the number of each amino acids among matched ones

if (@ARGV < 1) {die $usage;}

open(FA, "$ARGV[0]") || die $usage;

$seq = "";
while (<FA>) {
    if (/>\s*(.+)/) {
	if ($seq) {
	    @sym = split(//, $seq);
	    $ct = 0;
	    foreach $sym (@sym) {
		if ($sym eq "A") {
		    $nu{0} ++;
		    $ct ++;
		}
		if ($sym eq "B") {
		    $nu{1} ++;
		    $ct ++;
		}
		if ($sym eq "C") {
		    $nu{2} ++;
		    $ct ++;
		}
		if ($sym eq "D") {
		    $nu{3} ++;
		    $ct ++;
		}
		if ($sym eq "E") {
		    $nu{4} ++;
		    $ct ++;
		}
		if ($sym eq "F") {
		    $nu{5} ++;
		    $ct ++;
		}
		if ($sym eq "G") {
		    $nu{6} ++;
		    $ct ++;
		}	    
		if ($sym eq "H") {
		    $nu{7} ++;
		    $ct ++;
		}
		if ($sym eq "I") {
		    $nu{8} ++;
		    $ct ++;
		}
		if ($sym eq "K") {
		    $nu{9} ++;
		    $ct ++;
		}
		if ($sym eq "L") {
		    $nu{10} ++;
		    $ct ++;
		}
		if ($sym eq "M") {
		    $nu{11} ++;
		    $ct ++;
		}
		if ($sym eq "N") {
		    $nu{12} ++;
		    $ct ++;
		}
		if ($sym eq "P") {
		    $nu{13} ++;
		    $ct ++;
		}	    
		if ($sym eq "Q") {
		    $nu{14} ++;
		    $ct ++;
		}
		if ($sym eq "R") {
		    $nu{15} ++;
		    $ct ++;
		}
		if ($sym eq "S") {
		    $nu{16} ++;
		    $ct ++;
		}
		if ($sym eq "T") {
		    $nu{17} ++;
		    $ct ++;
		}
		if ($sym eq "V") {
		    $nu{18} ++;
		    $ct ++;
		}	    
		if ($sym eq "W") {
		    $nu{19} ++;
		    $ct ++;
		}
		if ($sym eq "Y") {
		    $nu{20} ++;
		    $ct ++;
		}
		if ($sym eq "Z") {
		    $nu{21} ++;
		    $ct ++;
		}
	}
	}
	if ($ct) {
	    printf "%05d ",$ct;
	}
	    &report;
	    $i = 0;
	    while ($i < 22) {
		$nu{$i} = 0;
		$i ++;
	    }
	    printf "\n>%s ", $1;
	    $seq = "";
    }
     else {
	chomp;
	$seq .= $_;
    }
}
close FA;

@sym = split(//, $seq);
$ct = 0;
foreach $sym (@sym) {
		if ($sym eq "A") {
		    $nu{0} ++;
		    $ct ++;
		}
		if ($sym eq "B") {
		    $nu{1} ++;
		    $ct ++;
		}
		if ($sym eq "C") {
		    $nu{2} ++;
		    $ct ++;
		}
		if ($sym eq "D") {
		    $nu{3} ++;
		    $ct ++;
		}
		if ($sym eq "E") {
		    $nu{4} ++;
		    $ct ++;
		}
		if ($sym eq "F") {
		    $nu{5} ++;
		    $ct ++;
		}
		if ($sym eq "G") {
		    $nu{6} ++;
		    $ct ++;
		}	    
		if ($sym eq "H") {
		    $nu{7} ++;
		    $ct ++;
		}
		if ($sym eq "I") {
		    $nu{8} ++;
		    $ct ++;
		}

		if ($sym eq "K") {
		    $nu{9} ++;
		    $ct ++;
		}
		if ($sym eq "L") {
		    $nu{10} ++;
		    $ct ++;
		}
		if ($sym eq "M") {
		    $nu{11} ++;
		    $ct ++;
		}
		if ($sym eq "N") {
		    $nu{12} ++;
		    $ct ++;
		}
		if ($sym eq "P") {
		    $nu{13} ++;
		    $ct ++;
		}	    
		if ($sym eq "Q") {
		    $nu{14} ++;
		    $ct ++;
		}
		if ($sym eq "R") {
		    $nu{15} ++;
		    $ct ++;
		}
		if ($sym eq "S") {
		    $nu{16} ++;
		    $ct ++;
		}
		if ($sym eq "T") {
		    $nu{17} ++;
		    $ct ++;
		}
		if ($sym eq "V") {
		    $nu{18} ++;
		    $ct ++;
		}	    
		if ($sym eq "W") {
		    $nu{19} ++;
		    $ct ++;
		}
		if ($sym eq "Y") {
		    $nu{20} ++;
		    $ct ++;
		}
		if ($sym eq "Z") {
		    $nu{21} ++;
		    $ct ++;
		}
	}
	
	    printf "%05d ",$ct;
	    &report;
	    $i = 0;
	    while ($i < 22) {
		$nu{$i} = 0;
		$i ++;
	    }
print "\n";

sub report {

    foreach $key (sort {$nu{$b} <=> $nu{$a}} (keys %nu)) {
	printf "%02d ",$nu{$key};
    }
}
