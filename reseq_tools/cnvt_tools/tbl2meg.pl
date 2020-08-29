use strict; 
use warnings; 
use SNP_tbl; 

!@ARGV and die "perl $0 in.snp\n"; 

my $inF = shift; 
my $st = SNP_tbl->new(filename=>$inF); 
$st->readTbl(); 
$st->tbl2meg(ofile=>"$inF.meg"); 


