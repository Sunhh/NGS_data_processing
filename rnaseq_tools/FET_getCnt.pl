#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"digits:i", # 0
	"addPref:s", # ''
	"outTotalCnt:s", 
#	"want:s", # 'F1' , or 'PP'
); 
$opts{'digits'} //= 0; 
$opts{'addPref'} //= ''; 
my $fh_totalCnt; 
if (defined $opts{'outTotalCnt'}) {
	$fh_totalCnt = &openFH( $opts{'outTotalCnt'}, '>' ); 
}

-t and !@ARGV and die "\nperl $0 TotalCnt_table.tbl.report_ralign_stat.fix P1g_perGene.sense.cnt.2_noSum_keep > P1g_perGene.sense.cnt.2_noSum_keep.forFET\n\n"; 

my $file_totalCnt = shift; 
my $file_cnt = shift; 

# Read in total mapped reads count. 
my %totalCnt; 
{
	my $fh = &openFH($file_totalCnt, '<'); 
	while (&wantLineC($fh)) {
		chomp; 
		my ($sample, $na1, $mapTtl) = &splitL("\t", $_); 
		$totalCnt{$sample} = $mapTtl; 
	}
	close($fh); 
}

# Generate data for FET analysis. 
my %h; 
my $fh_cnt = &openFH( $file_cnt, '<' ); 
CNT_LINE: 
while (<$fh_cnt>) {
	chomp; 
	my @ta = &splitL("\t", $_); 
	if ($. == 1) {
		for (my $i=1; $i<@ta; $i++) {
			my $sample = $ta[$i]; 
			$sample =~ m!^(P[13]g)_(allMisR0p\d+|P1spec|P3spec|cmmn|total)_([SM])(FL|FR|LV|RT|SD|ST)(F1|P1|P3)_(rep\d+)$! or die "Bad sample ID: [$sample]\n"; 
			my ($pG, $origin, $goodTag, $tissue, $indiv, $rep) 
			=  ($1,  $2,      $3,       $4,      $5,     $6); 
			$goodTag eq 'M' and next; 
			$goodTag eq 'S' or die "Bad goodTag [$goodTag] in [$sample]\n"; 
			$origin =~ m/^(P1spec|P3spec)$/ or next; 
			$indiv eq 'F1' or ($indiv eq 'P1' and $origin eq 'P1spec') or ($indiv eq 'P3' and $origin eq 'P3spec') or next; 
			$h{'tissue2idx'}{$tissue}{$indiv}{$origin}{$rep} = [$i, $sample, $totalCnt{$sample}]; 
		}
		## Clean tissue and samples. 
		my %new_tissue2idx; 
		for my $tissue (sort keys %{$h{'tissue2idx'}}) {
			# Good tissue must have related repID in P1/P3 and data in all individuals. 
			( defined $h{'tissue2idx'}{$tissue}{'F1'}{'P1spec'} and defined $h{'tissue2idx'}{$tissue}{'F1'}{'P3spec'} and defined $h{'tissue2idx'}{$tissue}{'P1'}{'P1spec'} and defined $h{'tissue2idx'}{$tissue}{'P3'}{'P3spec'} ) or next; 
			for my $rep (sort keys %{$h{'tissue2idx'}{$tissue}{'P1'}{'P1spec'}}) {
				defined $h{'tissue2idx'}{$tissue}{'P3'}{'P3spec'}{$rep} or next; 
				$new_tissue2idx{$tissue}{'P1'}{'P1spec'}{$rep} = [ @{ $h{'tissue2idx'}{$tissue}{'P1'}{'P1spec'}{$rep} } ]; 
				$new_tissue2idx{$tissue}{'P3'}{'P3spec'}{$rep} = [ @{ $h{'tissue2idx'}{$tissue}{'P3'}{'P3spec'}{$rep} } ]; 
			}
			defined $new_tissue2idx{$tissue}{'P1'}{'P1spec'} or next; 
			for my $rep (sort keys %{$new_tissue2idx{$tissue}{'P1'}{'P1spec'}}) {
				push( @{$h{'FETidx_P1P1s'}{$tissue}}, $new_tissue2idx{$tissue}{'P1'}{'P1spec'}{$rep}[0] ); 
				$h{'FETcnt_P1P1s'}{$tissue} += $new_tissue2idx{$tissue}{'P1'}{'P1spec'}{$rep}[2]; 
				push( @{$h{'FETidx_P3P3s'}{$tissue}}, $new_tissue2idx{$tissue}{'P3'}{'P3spec'}{$rep}[0] ); 
				$h{'FETcnt_P3P3s'}{$tissue} += $new_tissue2idx{$tissue}{'P3'}{'P3spec'}{$rep}[2]; 
			}
			@{$h{'FETidx_F1P1s'}{$tissue}} = map { $h{'tissue2idx'}{$tissue}{'F1'}{'P1spec'}{$_}[0] } sort keys %{$h{'tissue2idx'}{$tissue}{'F1'}{'P1spec'}}; 
			map { $h{'FETcnt_F1P1s'}{$tissue} += $h{'tissue2idx'}{$tissue}{'F1'}{'P1spec'}{$_}[2] } sort keys %{$h{'tissue2idx'}{$tissue}{'F1'}{'P1spec'}}; 
			@{$h{'FETidx_F1P3s'}{$tissue}} = map { $h{'tissue2idx'}{$tissue}{'F1'}{'P3spec'}{$_}[0] } sort keys %{$h{'tissue2idx'}{$tissue}{'F1'}{'P3spec'}}; 
			map { $h{'FETcnt_F1P3s'}{$tissue} += $h{'tissue2idx'}{$tissue}{'F1'}{'P3spec'}{$_}[2] } sort keys %{$h{'tissue2idx'}{$tissue}{'F1'}{'P3spec'}}; 
		}
		$h{'good_tissues'} = [ sort keys %{$h{'FETidx_P1P1s'}} ]; 
		for my $tissue (@{$h{'good_tissues'}}) {
			if ( $h{'FETcnt_P1P1s'}{$tissue} <= $h{'FETcnt_P3P3s'}{$tissue} ) {
				$h{'small_PP'}{$tissue} = 'P1P1s'; 
				$h{'high_PP'}{$tissue} = 'P3P3s'; 
				$h{'adj2small'}{$tissue} = $h{'FETcnt_P1P1s'}{$tissue}/$h{'FETcnt_P3P3s'}{$tissue}; 
			} else {
				$h{'small_PP'}{$tissue} = 'P3P3s'; 
				$h{'high_PP'}{$tissue} = 'P1P1s'; 
				$h{'adj2small'}{$tissue} = $h{'FETcnt_P3P3s'}{$tissue}/$h{'FETcnt_P1P1s'}{$tissue}; 
			}
		}
		my @out; 
		push(@out, "Gene_ID"); 
		for my $tissue (@{$h{'good_tissues'}}) {
			for my $fetS (qw/F1P1s F1P3s P1P1s P3P3s/) {
				my $hID = $opts{'addPref'} . "${tissue}_${fetS}_x" . scalar( @{$h{"FETidx_${fetS}"}{$tissue}} ); 
				push(@out, $hID); 
				defined $opts{'outTotalCnt'} and print {$fh_totalCnt} join("\t", $hID, 'Self_Count', $h{"FETcnt_${fetS}"}{$tissue})."\n"; 
			}
		}
		print STDOUT join("\t", @out)."\n"; 
		next CNT_LINE; 
	}
	my %cnt; 
	my @out; 
	push(@out, $ta[0]); 
	for my $tissue (@{$h{'good_tissues'}}) {
		for my $fetS (qw/F1P1s F1P3s P1P1s P3P3s/) {
			for my $idx (@{$h{"FETidx_${fetS}"}{$tissue}}) {
				$cnt{$tissue}{$fetS} += $ta[$idx]; 
			}
		}
		my $high_fetS = $h{'high_PP'}{$tissue}; 
		$cnt{$tissue}{$high_fetS} = sprintf( "%0.$opts{digits}f", $cnt{$tissue}{$high_fetS} * $h{'adj2small'}{$tissue} ); 
		push(@out, @{$cnt{$tissue}}{qw/F1P1s F1P3s P1P1s P3P3s/} ); 
	}
	print STDOUT join("\t", @out)."\n"; 
}
close($fh_cnt); 

