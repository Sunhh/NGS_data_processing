#!/usr/bin/perl
use strict; 
use warnings; 

!@ARGV and die "perl $0 augusts_test1.stdout augusts_test2.stdout > all_accuracy.tbl\n"; 

for my $fn (@ARGV) {
	my $v1 = &accuracy_calculator($fn); 
	print STDOUT "$fn\t$v1\n"; 
}

# calculate the result of testing AUGUSTUS on genbank files in a single number
# Copied from braker.pl
sub accuracy_calculator{
    my $aug_out=shift;
    my ($nu_sen, $nu_sp, $ex_sen, $ex_sp, $gen_sen, $gen_sp);
    open(AUGOUT, "$aug_out") or die ("Could not open $aug_out!\n");
    while(<AUGOUT>){
        if(/^nucleotide level\s*\|\s*(\S+)\s*\|\s*(\S+)/){
            $nu_sen=$1;
            $nu_sp=$2;
        }
        if(/^exon level\s*\|.*\|.*\|.*\|.*\|.*\|\s*(\S+)\s*\|\s*(\S+)/){
            $ex_sen=$1;
            $ex_sp=$2;
        }
        if(/^gene level\s*\|.*\|.*\|.*\|.*\|.*\|\s*(\S+)\s*\|\s*(\S+)/){
            $gen_sen=$1;
            $gen_sp=$2;
        }
    }
    my $target=(3*$nu_sen+2*$nu_sp+4*$ex_sen+3*$ex_sp+2*$gen_sen+1*$gen_sp)/15;
    return $target;
}

