#!/usr/bin/perl
# Trim */X in protein from tail. 
# Remove definition. 
# 50 characters per line. 
# Remove N in dna and remove X* in protein. 
use strict; 
use warnings; 
use fastaSunhh; 
use fileSunhh; 
my $fs_obj = fastaSunhh->new(); 

@ARGV or die "perl $0 pref\n"; 

my $fn_prot = "$ARGV[0].prot.fa"; 
my $fn_cds  = "$ARGV[0].cds.fa"; 

my $oh_prot = &openFH( "../$fn_prot", '>' ); 
my $oh_cds  = &openFH( "../$fn_cds" , '>' ); 

-e $fn_prot or $fn_prot = "$fn_prot.gz"; 
-e $fn_cds  or $fn_cds  = "$fn_cds.gz" ; 

my %seq_prot = %{ $fs_obj->save_seq_to_hash( 'faFile' => $fn_prot ) }; 
my %seq_cds  = %{ $fs_obj->save_seq_to_hash( 'faFile' => $fn_cds  ) }; 

scalar(keys %seq_prot) == scalar(keys %seq_cds) or die "diff\n"; 

print STDOUT join("\t", qw/ID cds_len prot_len Org.ID/)."\n"; 
for my $tk (sort { $seq_cds{$a}{'Order'} <=> $seq_cds{$b}{'Order'} } keys %seq_cds) {
	defined $seq_prot{$tk} or die "tk=$tk\n"; 
	$seq_prot{$tk}{'seq'} =~ s!\s!!g; 
	$seq_prot{$tk}{'seq'} =~ s![\*X]+$!!i; 
	$seq_cds{$tk}{'seq'}  =~ s!\s!!g; 
	$seq_cds{$tk}{'seq'}  =~ s![nN]!!i; 

	my $l_prot = length( $seq_prot{$tk}{'seq'} ); 
	my $l_cds  = length( $seq_cds{$tk}{'seq'} ); 
	unless ( $l_prot * 3 + 3 == $l_cds ) {
		print STDOUT "$tk\t$l_cds\t$l_prot\t$ARGV[0].$tk\n"; 
	}
	$l_prot > 0 or next; 
	$l_cds  > 0 or next; 

	$seq_prot{$tk}{'seq'} =~ s!(.{50})!$1\n!g; 
	chomp($seq_prot{$tk}{'seq'}); 
	$seq_cds{$tk}{'seq'}  =~ s!(.{50})!$1\n!g; 
	chomp($seq_cds{$tk}{'seq'} ); 
	print {$oh_prot} ">$ARGV[0].$tk\n$seq_prot{$tk}{'seq'}\n"; 
	print {$oh_cds}  ">$ARGV[0].$tk\n$seq_cds{$tk}{'seq'}\n"; 
}
close($oh_prot); 
close($oh_cds); 


