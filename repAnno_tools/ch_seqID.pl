#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and &stopErr("perl $0 in.out99.fa in.out99.dgt.tab 1>out99.named.fa 2>out99.new.tab\n"); 

my $inFaF = shift; 
my $inTabF = shift; 

open T,'<',"$inTabF" or die; 
# eleID   eleS    eleE    seqID   LTR1_S  LTR1_E  LTR2_S  LTR2_E  Inner_S Inner_E PBS_S   PBS_E   Strand
# eleID   eleS    eleE    Str     seqID   LTR1_S  LTR1_E  LTR2_S  LTR2_E  Inner_S Inner_E PBS_S   PBS_E   PPT_S   PPT_E
# RR1     98334   105001  ?       seq10   98339   99198   104137  104996  99199   104136  -1      -1      -1      -1
# RR2     785204  788754  +       seq10   785209  785369  788589  788749  785370  788588  785373  785384  -1      -1
my %region_to_ID; 
my (%RR2line, @RR_ids); 
while (<T>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my $tk1 = "$ta[4]_$ta[1]_$ta[2]"; # Element region. (with TSD)
	my $tk2 = "$ta[4]_$ta[5]_$ta[8]"; # LTR region. (without TSD)
	my $tk3 = "$ta[4]_$ta[9]_$ta[10]"; # Internal region. (without ltr region)
	defined $region_to_ID{$tk1} and die "$region_to_ID{$tk1}\n"; 
	defined $region_to_ID{$tk2} and die "$region_to_ID{$tk2}\n"; 
	defined $region_to_ID{$tk3} and die "$region_to_ID{$tk3}\n"; 
	$region_to_ID{$tk1} = "$ta[0]_${tk1}_ELE"; 
	$region_to_ID{$tk2} = "$ta[0]_${tk2}_LTR"; 
	$region_to_ID{$tk3} = "$ta[0]_${tk3}_INN"; 
	push(@RR_ids, $ta[0]); 
	$RR2line{$ta[0]} = [@ta]; 
}
close T; 

my (@in_id, @in_seq); 
my $is_oSeq=0; 
open F,'<',"$inFaF" or die; 
while (<F>) {
	if (m/^>/) {
		m/^>(\S+) \(dbseq\-nr (\d+)\) \[(\d+),(\d+)\]$/ or die "$_\n"; 
		my $seqID = "seq$2";
		my $scfID = "$1";
		my ($eleS, $eleE) = ($3,$4); 
		my $tk = "${seqID}_${eleS}_${eleE}"; 
		if ( defined $region_to_ID{$tk} ) {
			$is_oSeq = 1; 
			my $new_ID = "$region_to_ID{$tk}_$scfID"; 
			print STDOUT ">$new_ID\n"; 
			$new_ID =~ m/^(RR\d+)/ or die "$new_ID\n"; 
			my $rrID = $1; 
			push(@{$RR2line{$rrID}}, $scfID); 
		} else {
			$is_oSeq = 0; 
			# chomp($_); 
			# &tsmsg("[Rec] [Wrn] Record found in gff: tk=$tk $_\n"); 
		}
	}else{
		$is_oSeq == 1 and print STDOUT "$_"; 
	}
}
close F; 

for (@RR_ids) {
	$_ eq 'eleID' and push(@{$RR2line{$_}}, "scfID"); 
	print STDERR join("\t", @{$RR2line{$_}})."\n"; 
}
