#!/usr/bin/perl-w
use strict;
use Statistics::Descriptive; 

@ARGV == 1 or die "perl $0 struct_result_list\n"; 

#############Befor  Start  , open the files ####################
open (Res,"$ARGV[0]") || die "input file can't open $!";
my %OUT ;
################ Do what you want to do #######################
while($_=<Res>)
{
    chomp ;
    open (IA,"<$_") or die "$!" ;
# $/="Inferred clusters";
#Label (%Miss) :  Inferred clusters
# <IA> ;
# $/="\n";
# $_=<IA>;
my $K=-1; 
my $e='';
my $m='';
my $v=''; 
	while($_=<IA>) 
	{ 
        # 9     T_09    (3)   :  0.057 0.943
        #10     T_10    (0)   :  0.049 0.951
        #11     T_11    (0)   :  0.022 0.978
        chomp ; 
	if ($_ =~ m/^\s+(\d+)\s+populations assumed/) {
		$K = $1; 
	} elsif ( $_ =~ m/^Estimated Ln Prob of Data\s*=\s*([+-\d.]+)/ ) {
		push(@{$OUT{$K}{ELPD}}, $1); 
		$OUT{$K}{ELPD_SUM} += $1; 
		$OUT{$K}{ELPD_CNT} ++; 
		$e = $1; $m=$v=''; 
	} elsif ( $_ =~ m/^Mean value of ln likelihood\s*=\s*([+-\d.]+)/ ) {
		push(@{$OUT{$K}{MLPD}}, $1); 
		$OUT{$K}{MLPD_SUM} += $1; 
		$OUT{$K}{MLPD_CNT} ++; 
		$m=$1; $v=''; 
	} elsif ( $_ =~ m/^Variance of ln likelihood\s*=\s*([+-\d.]+)/ ) {
		$v=$1; 
		my $l = $m-$v/2; 
		push(@{$OUT{$K}{LK}}, $l); 
		$OUT{$K}{LK_SUM} += $l; 
		$OUT{$K}{LK_CNT} ++; 
		$K=-1; 
		$e=$m=$v=''; 
		last; 
	}
	}
close IA;
}
close Res ;
my @ak = sort {$a<=>$b} keys %OUT; 
for (my $i=0; $i<@ak; $i++) {
	my $K = $ak[$i]; 
	$OUT{$K}{'deltK'} = 'NA'; 
	if ($i >= 1 and $i < $#ak ) {
		my $k0 = $K-1; 
		my $k1 = $K; 
		my $k2 = $K+1; 
		my $stat2 = Statistics::Descriptive::Full->new(); 
		$stat2->add_data(@{$OUT{$k1}{LK}});
		my $sd_LK = $stat2->standard_deviation(); 
		my $stat1 = Statistics::Descriptive::Full->new(); 
		
		for (my $j=0; $j<@{$OUT{$k1}{LK}}; $j++) {
			my $v1 = abs($OUT{$k2}{LK}[$j]-2*$OUT{$k1}{LK}[$j]+$OUT{$k0}{LK}[$j]); 
			$stat1->add_data( $v1 ); 
		}
		my $v1 = $stat1->mean(); 
		$OUT{$K}{'deltK'} = $v1/$sd_LK; 
	}
	$OUT{$K}{'ELPD_AVG'} = (defined $OUT{$K}{ELPD_CNT}) ? $OUT{$K}{ELPD_SUM}/$OUT{$K}{ELPD_CNT} : 'NA' ; 
	$OUT{$K}{'MLPD_AVG'} = (defined $OUT{$K}{MLPD_CNT}) ? $OUT{$K}{MLPD_SUM}/$OUT{$K}{ELPD_CNT} : 'NA' ; 
	$OUT{$K}{'LK_AVG'}   = (defined $OUT{$K}{LK_CNT})   ? $OUT{$K}{LK_SUM}/$OUT{$K}{LK_CNT}     : 'NA' ; 
	$OUT{$K}{'K_p1'} //= 'NA'; 
	$OUT{$K}{'K_p2'} //= 'NA'; 
	$i >= 1 and $OUT{$K}{'K_p1'} = $OUT{$K}{'LK_AVG'} - $OUT{$K-1}{'LK_AVG'}; 
	$i >= 2 and $i < $#ak and $OUT{$K-1}{'K_p2'} = abs( $OUT{$K}{'K_p1'} - $OUT{$K-1}{'K_p1'} ); 

	$OUT{$K}{'runN'} = join(';', $OUT{$K}{'LK_CNT'}, $OUT{$K}{'ELPD_CNT'}, $OUT{$K}{'MLPD_CNT'}); 

	my $K_p1 = 'NA'; 
	my $K_p2 = 'NA'; 
	$i >= 1 and $K_p1 = $OUT{$K}{'LK_AVG'} - $OUT{$K-1}{'LK_AVG'}; 
	$OUT{$K}{'K_p1'} = $K_p1; 
	$i >= 2 and $K_p2 = abs(  ); 

}

print STDOUT join("\t", qw/K_NUM RUN_NUM deltK Avg_LK K_p1 K_p2 Avg_EstLPD Avg_AvgLPD/)."\n"; 
for (my $i=0; $i<@ak; $i++) {
	my $K = $ak[$i]; 
	my %oo = %{$OUT{$K}}; 
	print STDOUT join("\t", $K, @oo{qw/runN deltK LK_AVG K_p1 K_p2 ELPD_AVG MLPD_AVG/})."\n"; 
}



