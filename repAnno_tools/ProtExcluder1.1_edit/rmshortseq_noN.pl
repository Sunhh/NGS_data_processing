#! /usr/bin/perl
 
$usage = "rmshortseq.pl  stfile fastafile minmusize\n";

# to delete sequences short than given size

if (@ARGV < 3) {die "$usage";}

open(LENGTH, "$ARGV[0]") || die $usage;

while (<LENGTH>){
    if (/^\s*(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s/) {
	$nonN = $2 + $3 + $4 + $5;
	if ($nonN >= $ARGV[2]) {
	    $name = $1; 
     $len{$name} = $nonN; }
 }
}
close LENGTH;

open(FASTA, "$ARGV[1]") || die $usage;
while (<FASTA>) {
  
if   (/^>(\S+)s*/)  {
   
if (&comparison){$take = 1;}
    else {$take = 0;}
}
              if ($take){
		  print;
}
}
close FASTA;
                                                                                
sub comparison  {
    foreach $key (keys %len){
  if ($key eq $1)
  { return 1;}
}
}

