#!/usr/bin/perl
use strict; 
use warnings; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 
use LogInforSunhh; 

!@ARGV and die "perl $0 in.fas > out.meg\n"; 

my $file = shift; 
my %s2h = %{ $fs_obj->save_seq_to_hash( 'faFile'=>$file, 'has_head'=>1 ) }; 
print STDOUT <<HEADER; 
#mega
!Title tbl2meg;

HEADER

for my $tk (sort { $s2h{$a}{'Order'} <=> $s2h{$b}{'Order'} } keys %s2h) {
	$s2h{$tk}{'seq'} =~ s/\s//gs; 
	$s2h{$tk}{'seq'} =~ s/(.{60})/$1\n/g; chomp( $s2h{$tk}{'seq'} ); 
	print STDOUT "#$s2h{$tk}{'key'}\n$s2h{$tk}{'seq'}\n"; 
}

