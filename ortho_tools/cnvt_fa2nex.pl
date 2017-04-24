#!/usr/bin/perl
use strict; 
use warnings; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"datatype:s", # DNA
	"coding123!", # 
); 
$opts{'datatype'} //= 'DNA'; # Could be 'PROTEIN'; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 

my $help_txt = <<HH; 
perl $0 in_aln.fa > in_aln.fa.nex

-datatype     [DNA] protein; 
-coding123    [Boolean]
HH

!@ARGV and die $help_txt; 

my $fas = shift; 
my %sh = %{ $fs_obj->save_seq_to_hash('faFile'=>$fas) }; 

my ($ntax, $nchar); 

$ntax = scalar(keys %sh); 
for my $k1 (sort {$sh{$a}{'Order'} <=> $sh{$b}{'Order'}} keys %sh) {
	$sh{$k1}{'seq'} =~ s!\s!!g; 
	# $sh{$k1}{'seq'} =~ tr/N/?/; 
	$nchar //= length($sh{$k1}{'seq'}); 
}

print STDOUT <<HH;
#NEXUS

BEGIN DATA;
	DIMENSIONS NTAX=$ntax NCHAR=$nchar;
	FORMAT DATATYPE=$opts{'datatype'} GAP=-;

	MATRIX
HH

for my $k1 (sort {$sh{$a}{'Order'} <=> $sh{$b}{'Order'}} keys %sh) {
	print STDOUT "\t$k1\t$sh{$k1}{'seq'}\n"; 
}
print STDOUT ";\nEND;\n"; 

if ($opts{'coding123'}) {
	print <<HH; 
Begin assumptions;
	charset 1stpos = 1-$nchar\\3;
	charset 2ndpos = 2-$nchar\\3;
	charset 3rdpos = 3-$nchar\\3;
END;
HH
}
