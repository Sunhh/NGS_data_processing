#!/usr/bin/perl
use strict; 
use warnings; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 

!@ARGV and die "perl $0 qry_scf.fa ref_scf.fa V1ctg_to_hic.bn6.p100.u.ctgCov\n"; 

my $qsf = shift; 
my $tsf = shift; 

# ScfID   LenInCtg        CovLenInScf     CovPercInScf    CovByScf        CovLenSpan2     CovPercSpan2    CovS2   CovE2   CovS1   CovE1   LenInScf1       strands
# Cmo_Scf00001    11030317        11019172        99.90   HiC_scaffold_13 11019172        97.87   712309  11971090        1       11258782        11258782        +:11019172;-:
# Cmo_Scf00002    7545762 7541080 99.94   HiC_scaffold_16 7541080 98.05   14165614        21856752        1       7691139 7691139 +:7541080;-:0
# Cmo_Scf00003    7093495 7072078 99.70   HiC_scaffold_4  7072078 98.05   6942135 14155198        1       7213064 7213064 +:0;-:7072078
# Cmo_Scf00004    7060937 7058308 99.96   HiC_scaffold_10 7058308 98.29   1       7180833 1       7180833 7180833 +:0;-:7058308
# Cmo_Scf00005    6859295 6859295 100.00  HiC_scaffold_20 6859295 97.82   1       7012505 1       7012505 7012505 +:0;-:6859295

my %qseq = %{ $fs_obj->save_seq_to_hash('faFile' => $qsf) }; 
for (keys %qseq) { $qseq{$_}{'seq'} =~ s!\s!!g; $qseq{$_}{'seq'} = uc($qseq{$_}{'seq'}); $qseq{$_}{'len'} = length($qseq{$_}{'seq'});  }
my %tseq = %{ $fs_obj->save_seq_to_hash('faFile' => $tsf) }; 
for (keys %tseq) { $tseq{$_}{'seq'} =~ s!\s!!g; $tseq{$_}{'seq'} = uc($tseq{$_}{'seq'}); $tseq{$_}{'len'} = length($tseq{$_}{'seq'});  }

my %qhas; 

while (<>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	if (m!^\s*#! or $ta[0] eq 'ScfID') {
		print "$_\n"; 
		next; 
	}
	my ($qid,$qs,$qe, $tid, $ts, $te, $str1) = @ta[0,9,10, 4,7,8,12]; 
	my $str = &getStr($str1); 
	$qhas{$qid} = 1; 
	if ($str eq "+") {
		if ($qs != 1) {
			my $qsubSeq = substr($qseq{$qid}{'seq'}, 0, $qs-1); 
			my $tsubSeq = substr($tseq{$tid}{'seq'}, $ts-1-($qs-1), $qs-1); 
			if ($qsubSeq eq $tsubSeq) {
				$ts = $ts-($qs-1); 
				$qs = 1; 
			}
		}
		if ($qe != $qseq{$qid}{'len'}) {
			my $qsubSeq = substr($qseq{$qid}{'seq'}, $qe, $qseq{$qid}{'len'}-$qe); 
			my $tsubSeq = substr($tseq{$tid}{'seq'}, $te, $qseq{$qid}{'len'}-$qe); 
			if ($qsubSeq eq $tsubSeq) {
				$te = $te+$qseq{$qid}{'len'}-$qe; 
				$qe = $qseq{$qid}{'len'}; 
			}
		}
	} elsif ($str eq "-") {
		my ($ori_ts, $ori_te) = ($ts, $te); 
		if ($qs != 1) {
			my $qsubSeq = substr($qseq{$qid}{'seq'}, 0, $qs-1); 
			my $tsubSeq = substr($tseq{$tid}{'seq'}, $ori_te, $qs-1); 
			$tsubSeq =~ tr/acgturykmbvdhACGTURYKMBVDHwWsSnN/tgcaayrmkvbhdTGCAAYRMKVBHDwWsSnN/;
			$tsubSeq = reverse($tsubSeq); 
			if ($qsubSeq eq $tsubSeq) {
				$te = $ori_te+($qs-1); 
				$qs = 1; 
			}else {
				# die "q=$qsubSeq\nt=$tsubSeq\n"; 
				; 
			}
		}
		if ($qe != $qseq{$qid}{'len'}) {
			my $qsubSeq = substr($qseq{$qid}{'seq'}, $qe, $qseq{$qid}{'len'}-$qe); 
			my $tsubSeq = substr($tseq{$tid}{'seq'}, $ori_ts-1-($qseq{$qid}{'len'}-$qe), $qseq{$qid}{'len'}-$qe); 
			$tsubSeq =~ tr/acgturykmbvdhACGTURYKMBVDHwWsSnN/tgcaayrmkvbhdTGCAAYRMKVBHDwWsSnN/;
			$tsubSeq = reverse($tsubSeq); 
			if ($qsubSeq eq $tsubSeq) {
				$ts = $ori_ts-($qseq{$qid}{'len'}-$qe); 
				$qe = $qseq{$qid}{'len'}; 
			}
		}
	} else {
		die "bad str: $str\n"; 
	}
	if ( $qs==1 and $qe == $qseq{$qid}{'len'} ) {
		my $tsubSeq = substr($tseq{$tid}{'seq'}, $ts-1, $te-$ts+1); 
		if ($str eq '-') {
			$tsubSeq =~ tr/acgturykmbvdhACGTURYKMBVDHwWsSnN/tgcaayrmkvbhdTGCAAYRMKVBHDwWsSnN/; 
			$tsubSeq = reverse($tsubSeq); 
		}
		if ( $qseq{$qid}{'seq'} eq $tsubSeq ) {
			# This is good.
			@ta[9,10,7,8] = ($qs,$qe, $ts,$te); 
			print STDOUT join("\t", @ta)."\n"; 
		} else {
			@ta[9,10,7,8] = ($qs,$qe, $ts,$te); 
			print STDERR join("\t", @ta)."\n"; 
		}
	} else {
		print STDERR join("\t", @ta)."\n"; 
	}
}

# For those not covered
for my $k1 (keys %qseq) {
	defined $qhas{$k1} and next; 
	my @idx;
	for my $k2 (keys %tseq) {
		my $ii = index($tseq{$k2}{'seq'}, $qseq{$k1}{'seq'}); 
		$ii > -1 and push(@idx, [$k2, $ii+1, $ii+$qseq{$k1}{'len'}, "+"]); 
		my $s1 = reverse($qseq{$k1}{'seq'}); 
		$s1 =~ tr/acgturykmbvdhACGTURYKMBVDHwWsSnN/tgcaayrmkvbhdTGCAAYRMKVBHDwWsSnN/; 
		my $i2 = index($tseq{$k2}{'seq'}, $s1); 
		$i2 > -1 and push(@idx, [$k2, $i2+1, $i2+$qseq{$k1}{'len'}, "-"]); 
	}
	if (scalar(@idx) == 1) {
		print STDOUT join("\t", $k1, "x", "x", "x", $idx[0][0], "x", "x", $idx[0][1], $idx[0][2], 1, $qseq{$k1}{'len'}, $qseq{$k1}{'len'}, $idx[0][3])."\n"; 
	} elsif (scalar(@idx) > 1) {
		for my $i3 (@idx) {
			print STDERR join("\t", $k1, "m", "m", "m", $i3->[0], "m", "m", $i3->[1], $i3->[2], 1, $qseq{$k1}{'len'}, $qseq{$k1}{'len'}, $i3->[3])."\n"; 
		}
	} else {
		print STDERR join("\t", $k1, qw/u u u u u u u u u u/, $qseq{$k1}{'len'}, "+:0;-:0")."\n"; 
	}
}

sub getStr {
	my @a1;
	for my $t1 (split(/;/, $_[0])) {
		push(@a1, [split(/:/, $t1)]); 
	}
	@a1 = sort { $b->[1]<=>$a->[1] } @a1; 
	return($a1[0][0]); 
}
