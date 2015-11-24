#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 

my $help_txt = <<HH; 

perl $0 in_pca.set01_all.pca.evec > in_pca.set01_all.pca.evec.tbl 

HH

@ARGV >= 1 or &LogInforSunhh::usage($help_txt); 

my @evec_arr = @{ &load_pca_evec( $ARGV[0] ) }; 
# my %ind_info = %{ &load_pca_indv( $ARGV[1] ) }; 

print "Indv\tGrp"; 
for ( my $i=1; $i<@{$evec_arr[0]}; $i++) {
	print "\tEV$i"; 
}
print "\n"; 

for my $tr1 ( @evec_arr[ 1 .. $#evec_arr ] ) {
	print "$tr1->[0]\t$tr1->[-1]"; 
	for ( my $i=1; $i+1<@$tr1; $i++ ) {
		print "\t$tr1->[$i]"; 
	}
	print "\n"; 
}
&tsmsg("[Msg] Done for $0\n"); 

sub load_pca_evec {
	# $_[0]: in_pca.set01_all.pca.evec 
	my $fh = &openFH($_[0], '<'); 
	my @back; 
	my $maxColN; 
	while (<$fh>) {
		chomp; 
		s/^\s+//;
		my @ta = split(/\s+/, $_); 
		$maxColN //= $#ta; 
		push(@back, [@ta]); 
	}
	close($fh); 
	return(\@back); 
}# load_pca_evec 

sub load_pca_indv {
	# $_[0] : in_pca.set01_all.ind.shrt 
	my $fh = &openFH($_[0], '<'); 
	my %back; 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		if ( defined $back{$ta[0]} ) {
			&tsmsg("[Wrn] repeated ID [$ta[0]]\n"); 
			next; 
		}
		$back{$ta[0]} = [@ta[1 .. $#ta]]; 
	}
	close($fh); 
	return (\%back); 
}# load_pca_indv 

