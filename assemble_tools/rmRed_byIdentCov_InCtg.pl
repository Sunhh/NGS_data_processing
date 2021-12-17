#!/usr/bin/perl
# 12/15/2021 : Fix the bug when there are two exactly the same contigs in the file, both of them will be removed. 
use strict; 
use warnings; 
use LogInforSunhh; 


!@ARGV and die "perl $0 X20.ctg.fa [min_identityRate [min_coverageRate]]\n"; 

my $fn = shift; 
my $min_identR = 0.99; 
my $min_covR   = 0.99; 
scalar(@ARGV) > 0 and $min_identR = shift; 
scalar(@ARGV) > 0 and $min_covR = shift; 

my $max_seqlen = 3e5; # 300kb is a good length maximum to remove redundancy. 
my $bn_wordsize = 100; 
my $bn_eval     = 1e-10; 

my @torm; 

&runCmd("deal_fasta.pl -attr key:len $fn > $fn.kl"); push(@torm, "$fn.kl"); 

open F,'<',"$fn.kl" or die; 
open O1,'>',"$fn.kl.short_list" or die; push(@torm, "$fn.kl.short_list"); 
while (<F>) {
	chomp; 
	my @ta=split(); 
	$ta[1] eq 'len' and next; 
	$ta[1] < $max_seqlen and print O1 "$_\n"; 
}
close O1; 
close F; 


&runCmd("deal_fasta.pl $fn -drawByList -drawList $fn.kl.short_list -drawLcol 0 -drawWhole > $fn.kl.short_list.fa"); push(@torm, "$fn.kl.short_list.fa"); 
&runCmd("makeblastdb -dbtype nucl -in $fn"); push(@torm, "$fn.n??"); # I don't want to remove short sequences in the database because there might be short primary contigs that covers shorter redundancy. 
# &runCmd("deal_fasta.pl -scf2ctg $fn.kl.short_list.c2g $fn.kl.short_list.fa"); 
&runCmd("blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand' -num_threads 30 -word_size $bn_wordsize -evalue $bn_eval -db $fn -query $fn.kl.short_list.fa -out $fn.kl.short_list.fa.toAll.bn6"); push(@torm, "$fn.kl.short_list.fa.toAll.bn6"); 
open F2,'<',"$fn.kl.short_list.fa.toAll.bn6" or die; 
open O2,'>',"$fn.kl.short_list.fa.toAll.bn6.red" or die; 
# push(@torm, "$fn.kl.short_list.fa.toAll.bn6.red"); 
my %u2; 
while (<F2>) {
	chomp; 
	my @ta=split; 
	$ta[0] eq $ta[1] and next; 
	$ta[2] >= 100*$min_identR or next; 
	($ta[7]-$ta[6]+1) >= $min_covR * $ta[12] or next; 
	defined $u2{$ta[0]} and next; 
	if ($ta[12] < $ta[13]) {
		$u2{$ta[0]} = 1; 
	} elsif ($ta[12] == $ta[13]) {
		defined $u2{$ta[1]} or $u2{$ta[0]} = 1; # I need some more modifications here to get an accurate result.
	} else {
		defined $u2{$ta[0]} or $u2{$ta[1]} = 1;
	}
	print O2 "$_\n"; 
}
close O2; 
close F2; 

# &runCmd("deal_table.pl -kSrch_drop -kSrch_idx $fn.kl.short_list.fa.toAll.bn6.red $fn.kl.short_list > $fn.kl.short_list.noRed"); 
&runCmd("deal_table.pl -kSrch_drop -kSrch_idx $fn.kl.short_list.fa.toAll.bn6.red $fn.kl > $fn.kl.noRed"); 
&runCmd("deal_fasta.pl -drawByList -drawLcol 0 -drawWhole -drawList $fn.kl.noRed $fn > $fn.noRed.fa"); 
&runCmd("deal_fasta.pl -N50 $fn.noRed.fa > $fn.noRed.fa.N50"); 

for my $a (@torm) {
	system "rm -f $a"; 
}


sub runCmd {
	&exeCmd_1cmd($_[0]) and &stopErr("[Err] Failed at CMD: $_\n"); 
}
