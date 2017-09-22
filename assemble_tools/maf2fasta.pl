#!/usr/bin/perl
use strict; 
use warnings; 
use ReadInAlnSunhh; 
use Getopt::Long; 

my %opts; 
GetOptions(\%opts, 
	"help!", 
	"label:s", 
	"out:s", 
); 

-t and !@ARGV and &usage();
$opts{help} and &usage();

sub usage {
	print STDOUT <<HELP; 
perl $0 in.maf
-help 
-label     [''] Not used. 
HELP
	exit 1; 
}

my $oFh = \*STDOUT;
if (defined $opts{out}) {
	my $tfh; 
	open $tfh, '>', "$opts{out}" or die;
	$oFh = $tfh; 
}

my @FH; 
!(-t) and push (@FH, \*STDIN); 
for (@ARGV) {
	my $tfh; 
	open $tfh, '<', "$_" or die; 
	push(@FH, $tfh); 
}

my @all_blks; 
for my $fh (@FH) {
	while ( my %rec1 = %{ readMAF($fh) } ) {
		chomp( $rec1{a}[0] ); 
		my @lines; 
		for my $tline ( @{$rec1{o}} ) {
			$tline =~ m/^s\s/ or next; 
			chomp($tline); 
			push(@lines, $tline); 
		}
		scalar(@lines) >= 2 or next; 
		push(@all_blks, [$rec1{a}[0]]); 
		for my $tline (@lines) {
			my $tr = splitMafSline($tline, 1); 
			push(@{$all_blks[-1]}, $tr); 
		}
	}
}

for my $tr1 (sort { $a->[1]{seqId} cmp $b->[1]{seqId} || $a->[1]{seqStrand} cmp $b->[1]{seqStrand} || $a->[1]{seqStart} <=> $b->[1]{seqStart} } @all_blks) {
# for my $tr1 ( @all_blks) {
	for (my $i=1; $i<@$tr1; $i++) {
		my %tr2 = %{$tr1->[$i]}; 
		print {$oFh} ">$tr2{seqId} $tr2{seqStart} $tr2{blkSize} $tr2{seqStrand} $tr2{seqLen}\n"; 
		$tr2{seqSeq} =~ s/(.{60})/$1\n/g; 
		chomp($tr2{seqSeq}); 
		print {$oFh} $tr2{seqSeq} . "\n"; 
	}
	print {$oFh} "=\n"; 
}



