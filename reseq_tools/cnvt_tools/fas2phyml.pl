#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	
); 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 

my $help_txt = <<HH; 

perl $0 in_aligned.fasta > out_phyml_interleaved.phyml

-help       


HH

!@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

&tsmsg("[Msg] Reading fasta.\n\n"); 
my $charN_per_line = 60; 

my %fas_hash = %{$fs_obj->save_seq_to_hash('faFile'=>$ARGV[0], 'has_head'=>1)}; 
my @seq_IDs = sort { $fas_hash{$a}{'Order'} <=> $fas_hash{$b}{'Order'} } keys %fas_hash; 
for (@seq_IDs) {
	$fas_hash{$_}{'seq'} =~ s/\s//g; 
	$fas_hash{$_}{'len'} = length($fas_hash{$_}{'seq'}); 
}
my $seq_num = scalar(@seq_IDs); 
my $seq_len = $fas_hash{ $seq_IDs[0] }{'len'}; 

&tsmsg("[Msg] Printint output.\n"); 
print STDOUT "$seq_num $seq_len\n"; 

for (my $i=0; $i<$seq_len; $i+=$charN_per_line) {
	if ($i==0) {
		for (@seq_IDs) {
			my $new_ID = substr($_, 0, 10); 
			$new_ID =~ s/[\s();:]/_/g; 
			$new_ID = sprintf("%-9s", $new_ID); 
			print STDOUT join(' ', $new_ID, @{ &arr_by_seg( substr($fas_hash{$_}{'seq'}, $i, $charN_per_line), 10 ) })."\n"; 
		}
		print STDOUT "\n"; 
		next; 
	}
	for (@seq_IDs) {
		print STDOUT join(' ', ' ' x 9, @{ &arr_by_seg( substr($fas_hash{$_}{'seq'}, $i, $charN_per_line), 10 ) })."\n"; 
	}
	print STDOUT "\n"; 
}

sub arr_by_seg {
	# $_[0] whole sequence 
	# $_[1] length per segment. 
	my @back; 
	for (my $i=0; $i<length($_[0]); $i+=$_[1]) {
		my $ss = substr($_[0], $i, $_[1]); 
		push(@back, $ss); 
	}
	return \@back; 
}




