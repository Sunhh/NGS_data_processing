#!/usr/bin/perl
# 2019-11-06 Not considering the e-value currently. 
use strict; 
use warnings; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
  "ncol_score:i", # 11
  "ncol_index:i", # 0
  "max_diffR:f",  # 0
  "help!", 
); 

$opts{'ncol_score'} //= 11; 
$opts{'ncol_index'} //= 0; 
$opts{'max_diffR'}  //= 0; 

my $htxt = <<HH; 
#######################################
# perl $0  pep.bp6 > pep.bp6.slct
#
# -max_diffR     [$opts{'max_diffR'}] [0-1); 
#
# -ncol_score    [$opts{'ncol_score'}]
# -ncol_index    [$opts{'ncol_index'}]
#
HH
$opts{'help'} and die "$htxt"; 
-t and !@ARGV and die "$htxt"; 

my %ss; 
my $n=0;
my $has_bs = 0;
while (<>) {
  chomp; 
  my @ta=split(/\t/, $_); 
  if ($ta[ $opts{'ncol_score'} ] =~ m!^bitscore$!i) {
    $has_bs == 0 or next;
    $has_bs = 1;
    print STDOUT "$_\n";
    next;
  }
  if (defined $ss{$ta[$opts{'ncol_index'}]}{'topScore'}) {
    $ss{$ta[$opts{'ncol_index'}]}{'topScore'} < $ta[ $opts{'ncol_score'} ] and $ss{$ta[$opts{'ncol_index'}]}{'topScore'} = $ta[ $opts{'ncol_score'} ]; 
  } else {
    $ss{$ta[$opts{'ncol_index'}]}{'topScore'} = $ta[ $opts{'ncol_score'} ]; 
    $n++; 
    $ss{ $ta[$opts{'ncol_index'}] }{'rank'} = $n; 
  }
  $ta[ $opts{'ncol_score'} ] >= (1-$opts{'max_diffR'}) * $ss{$ta[$opts{'ncol_index'}]}{'topScore'} or next; 
  push(@{$ss{$ta[$opts{'ncol_index'}]}{'ele'}}, [ $_, $ta[ $opts{'ncol_score'} ] ]); 
}
for my $i1 (sort { $ss{$a}{'rank'} <=> $ss{$b}{'rank'} } keys %ss) {
  my @a1; 
  for my $e1 (@{$ss{$i1}{'ele'}}) {
    $e1->[1] >= $ss{$i1}{'topScore'} * (1-$opts{'max_diffR'}) or next; 
    push(@a1, $e1); 
  }
  @a1 = sort { $b->[1] <=> $a->[1] } @a1; 
  for my $e1 (@a1) {
    print STDOUT "$e1->[0]\n"; 
  }
}

