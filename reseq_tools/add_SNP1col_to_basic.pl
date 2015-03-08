#!/usr/bin/perl -w 
# 2013-08-19 Edit to ignore the 3rd column information in WM97_refChr file. 
use strict; 

my $snpLisF = shift; 
my $refChrF = shift; 
!@ARGV and die "perl $0 basicSNPlist refChr 1col.1 1col.2 ...\n"; 

my @save_idx; 
my $local_time; 

$local_time = localtime(); 
warn "[Stat] $local_time Program $0 started.\n"; 

{
 my @snp_data; 
 &fast_readF($snpLisF, \@snp_data); 
 $local_time = scalar(localtime()); 
 warn "[Stat] $local_time SNP list [$snpLisF] read in.\n"; 
 my %h; 
 for (my $i=0; $i<@snp_data; $i++) {
  $snp_data[$i] eq '' and next; 
#  $snp_data[$i] =~ m/^(\S+\t\S+\t\S+)/ or die "i=$i, $snp_data[$i]\n"; # OLD 
  $snp_data[$i] =~ m/^(\S+\t\S+)/ or die "i=$i, $snp_data[$i]\n"; 
  my $k = $1; 
  defined $h{$k} and die "repeat $k\n"; 
  $h{$k} = $i; 
 }
 $local_time = scalar(localtime()); 
 warn "[Stat] $local_time SNP list [$snpLisF] parsed.\n"; 
 my @ref_data; 
 &fast_readF($refChrF, \@ref_data); 
 $local_time = scalar(localtime()); 
 warn "[Stat] $local_time RefChr [$refChrF] read in.\n"; 
 for (my $i=0; $i<@ref_data; $i++) {
	 $i % 100000 == 0 and warn "[Msg]$i lines OK. " . scalar(localtime()) . "\n"; 
  chomp($ref_data[$i]); 
	$ref_data[$i] =~ m/^(\S+\t\S+)/ or die "i=$i, $ref_data[$i]\n"; 
	my $k = $1; 
#  defined $h{$ref_data[$i]} and push( @save_idx, [ $i, $snp_data[$h{$ref_data[$i]}] ] ); # OLD 
  defined $h{$k} and push( @save_idx, [ $i, $snp_data[$h{$k}] ] ); 
 }
 $local_time = scalar(localtime()); 
 warn "[Stat] $local_time Positions in RefChr recorded.\n"; 
}

for my $fn (@ARGV) {
 my @add_data; 
 &fast_readF($fn, \@add_data); 
 for my $tr (@save_idx) {
  $tr->[1] .= "\t$add_data[$tr->[0]]"; 
 }
 $local_time = scalar(localtime()); 
 warn "[Stat] $local_time SNP file $fn added.\n"; 
}

$local_time = scalar(localtime());
warn "[Stat] $local_time All stored. Printing.\n"; 
for (@save_idx) {
 print STDOUT "$_->[1]\n"; 
}

$local_time = scalar(localtime());
warn "[Stat] $local_time Finish.\n"; 

# store file lines in the second value; 
sub fast_readF {
 # local $/ = undef(); 
 if ($_[0] =~ m/\.gz$/) {
  open F,'-|',"gzip -cd $_[0]" or die; 
 } else {
  open F,'<',"$_[0]" or die; 
 }
 my @data; 
 while (<F>) {
	$. % 10000000 == 1 and warn "[Msg] $. lines [$_[0]] [" . scalar(localtime()) . "]\n"; 
  chomp; 
  push(@data, $_); 
 }
 # my $text = <F>; 
 close F; 
 # @{$_[1]} = split(/\n/, $text, -1); 
 @{$_[1]} = @data; 
 return 0; 
}

