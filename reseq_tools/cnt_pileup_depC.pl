#!/usr/bin/perl
use strict; 
use warnings; 

use LogInforSunhh; 

for (@ARGV) {
	&exeCmd("awk '\$3 != \"N\" { print \$4 }' $_ | deal_table.pl -col_repCount 0 | deal_table.pl -col_sort 1 > $_.noN.depC"); 
}
