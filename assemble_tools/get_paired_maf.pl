#!/usr/bin/perl
use strict; 
use warnings; 
use ReadInAlnSunhh; 
use Getopt::Long;

my %opts; 
GetOptions(\%opts, 
	"help!", 
	"out:s", 
); 

-t and !@ARGV and &usage(); 
$opts{help} and &usage(); 

sub usage {
print STDOUT <<HELP; 
perl $0 in.maf
-help 
-out        [\\*STDOUT] Output file name. 
HELP
	exit 1; 
}

my $oFh = \*STDOUT; 
if (defined $opts{out}) {
	my $tfh; 
	open $tfh,'>',"$opts{out}" or die "Failed to open $opts{out}\n$!\n"; 
	$oFh = $tfh; 
}

my @FH; 
!(-t) and push(@FH, \*STDIN); 
for (@ARGV) {
	my $fh; 
	open $fh, '<', "$_" or die; 
	push(@FH, $fh); 
}

print {$oFh} "##maf version=1\n"; 
my @all_blks; 
for my $fh (@FH) {
	while ( my %rec1 = %{readMAF($fh)} ) {
		chomp( $rec1{a}[0] ); 
		my @lines; 
		for my $tline (@{$rec1{o}}) {
			$tline =~ m/^s\s/ or next; 
			chomp($tline); 
			push(@lines, $tline); 
		}
		scalar(@lines) >= 2 or next; 
		print {$oFh} $rec1{a}[0] . "\n"; 
		for my $tline (@lines) {
			print {$oFh} $tline . "\n"; 
		}
		print {$oFh} "\n"; 
	}
}


