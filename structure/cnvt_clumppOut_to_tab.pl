#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"clumpp_out:s", # ClumppIndFile.output 
	"new_order:s",  # A file with format: ID_in_clumpp \\t ID_text , this also defined the order of output table. 
	"help!", 
); 

my $help_txt = <<HH; 

perl $0 -clumpp_out ClumppIndFile.output      -new_order ordered_ID_list 

Format of ordered_ID_list : ID_in_clumpp \\t ID_text 

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
for my $pp (qw/clumpp_out new_order/) {
	defined $opts{$pp} or &LogInforSunhh::usage($help_txt); 
}

my %clumpp_data = %{ &load_clumpp_out( $opts{'clumpp_out'} ) }; 
my @order_list  = @{ &load_id_list( $opts{'new_order'} ) }; 

my $has_header = 0; 
for my $a1 ( @order_list ) {
	my ( $id_clumpp, $id_text ) = @$a1;  
	defined $clumpp_data{ $id_clumpp } or &stopErr("[Err] Undefined ID_in_clumpp [$id_clumpp]\n"); 

	if ($has_header == 0) {
		print STDOUT "IndID"; 
		for (my $i=1; $i-1 < @{$clumpp_data{$id_clumpp}[0]}; $i++) { 
			print STDOUT "\tC_$i"; 
		}
		print STDOUT "\n"; 
		$has_header = 1; 
	}
	print STDOUT join("\t", $id_text, @{$clumpp_data{$id_clumpp}[0]})."\n"; 
}


sub load_id_list {
	# $_[0] : 
	my $fh = &openFH($_[0], '<'); 
	my @a; 
	while (<$fh>) {
		m/^\s*(#|$)/ and next; 
		chomp; 
		my @ta = split(/\t/, $_); 
		push(@a, [$ta[0], $ta[1]]); # [ID_in_clumpp, ID_text]
	}
	close($fh); 
	return(\@a); 
}# load_id_list() 

sub load_clumpp_out {
	# $_[0] : ClumppIndFile.output 
	my %h; 
	my $fh = &openFH($_[0], '<'); 
	my $colN ; 
	while (<$fh>) {
		chomp; 
		s/^\s+//; 
		s/\s+$//; 
		my @ta = split(/\s+/, $_); 
		$ta[4] eq ':' or &stopErr("[Err] Unknown format in ClumppIndFile.output : $_\n"); 
		$colN //= $#ta; 
		defined $h{$ta[1]} and &stopErr("[Err] Repeated indID in ClumppIndFile.output : $_\n"); 
		$h{$ta[1]} = [ [@ta[5 .. $colN]], $ta[0] ]; # $h{ ID_in_clumpp } = [ [coefficient1, coefficient2, ..., coefficientK], GrpID ]
	}
	close($fh); 
	return(\%h); 
}# load_clumpp_out() 


