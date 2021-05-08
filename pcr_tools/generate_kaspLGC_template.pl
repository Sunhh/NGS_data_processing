#!/usr/bin/perl
use strict; 
use warnings; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 

!@ARGV and die "perl $0 ref.fasta in_table > out_table\n"; 
# Example of in_table: 
#chr2    78732   C       T
#chr2    139549  T       C
#chr2    278407  G       T
#chr2    384663  T       A
#

my $flen = 300; 

my $fafile = shift; 
my %ss = %{ $fs_obj->save_seq_to_hash('faFile'=>$fafile) }; 
for (keys %ss) { $ss{$_}{'seq'} =~ s!\s!!g; $ss{$_}{'seq'} = uc($ss{$_}{'seq'}); $ss{$_}{'len'}=length($ss{$_}{'seq'}); }

while (<>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	unless (defined $ss{$ta[0]}) {
		print join("\t", @ta, "NA")."\n"; 
		next; 
	}
	my ($l_s, $l_e) = ( $ta[1]-1-$flen+1, $ta[1]-1 ); 
	$l_s < 1 and $l_s = 1; 
	my ($r_s, $r_e) = ( $ta[1]+1, $ta[1]+1+$flen-1 ); 
	$r_e > $ss{$ta[0]}{'len'} and $r_e = $ss{$ta[0]}{'len'}; 
	my $lseq = substr( $ss{$ta[0]}{'seq'}, $l_s-1, $l_e-$l_s+1 ); 
	my $rseq = substr( $ss{$ta[0]}{'seq'}, $r_s-1, $r_e-$r_s+1 ); 
	$lseq =~ s!^(.[NU]+)!! and $l_s -= length($1); 
	$rseq =~ s!([NU]+.)$!! and $l_e -= length($1); 
	if ($l_e - $l_s +1 < 15 or $r_e - $r_s + 1 < 15) {
		print join("\t", @ta, "FAILED")."\n"; 
		next; 
	}
	my $oseq = $lseq . "[$ta[2]/$ta[3]]" . $rseq; 
	print join("\t", @ta, $oseq)."\n"; 
}


