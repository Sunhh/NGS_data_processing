#!/usr/bin/perl
use strict; 
use warnings; 
use Text::Balanced qw( extract_bracketed ); 
# http://www.perlmonks.org/?node_id=547596
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"branch_scale:f", # 
	"branch_round:i", # -1 
); 
$opts{'branch_scale'} //= 1; 
$opts{'branch_round'} //= -1; 

my $help_txt = <<HH; 

perl $0 old2new_ID in.nwk > new.nwk 

-branch_scale    [$opts{'branch_scale'}]
-branch_round    [$opts{'branch_round'}]

HH

@ARGV == 2 or &LogInforSunhh::usage($help_txt); 

my %id2new = %{ &load_id2new($ARGV[0]) }; 
open AA,'<',"$ARGV[1]" or die; 
while (<AA>) {
	my $eee = chomp; 
	my @posList; 
	pos($_) = 0; 
	while ($_ =~ m/\G(?:.*?[\(,])([^\s\(\),:]+)(?:\:)/gs) {
		push( @posList, [ $-[1]+1, $+[1], $1 ] ); 
		pos($_) = $+[1]; 
	}
	my $prevE = 0; 
	my @ele; 
	for my $ar (@posList) {
		my ($cS, $cE, $match) = @$ar; 
		if ($cS-1 > $prevE) {
			push( @ele, substr($_, $prevE, $cS-1-$prevE) ); 
		}
		my $toChg = substr($_, $cS-1, $cE-$cS+1); 
		defined $id2new{$toChg} and $toChg = $id2new{$toChg}; 
		push( @ele, $toChg ); 
		$prevE = $cE; 
	}
	if ( length($_) >= $prevE+1 ) {
		push( @ele, substr($_, $prevE, length($_)-$prevE) ); 
	}
	for (my $i=0; $i<@ele; $i++) {
		$ele[$i] = &fmtBranch($ele[$i]); 
	}
	print join('', @ele) . ( ($eee > 0) ? "\n" : "" ) ;
}
close AA; 

sub load_id2new {
	my $fh = &openFH($_[0], '<'); 
	my %back; 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		$back{$ta[0]} = $ta[1]; 
	}
	close($fh); 
	return(\%back); 
}

#while (<>) {
#	chomp; 
#	my $n = &rmStatNum($_); 
#	my $m = &fmtBranch($n); 
#	print "$m\n"; 
#}

# (((Sly:0.11746590,Vvi:0.07461103)0.7300:0.02307856,Ath:0.14512370)1.0000:0.03408003,(SpiOl:0.05289341,Bvu:0.04303894)1.0000:0.05605601);

sub rmStatNum {
	while ($_[0] =~ s/\)\s*[\d.]+/)/g) {
	# while ($_[0] =~ s/:[\d.]+//g or $_[0] =~ s/\)\s*[\d.]+/)/g) {
	}
	return $_[0]; 
}
sub fmtBranch {
	if ( $opts{'branch_round'} >= 0 ) {
		$_[0] =~ s/:(\d+\.\d+)/":" . sprintf("%.0f", $1*$opts{'branch_scale'})/eg
		#my $prev = $_[0]; 
		# while ($_[0] =~ s/:(\d+\.\d+)/":" . int($1*$opts{'branch_scale'})/eg) {
		#while ($_[0] =~ s/:(\d+\.\d+)/":" . sprintf("%.0f", $1*$opts{'branch_scale'})/eg) {
		#	$prev eq $_[0] and last; 
		#	$prev = $_[0]; 
		#}
	} else {
		$_[0] =~ s/:(\d+\.\d+)/":" . ($1*$opts{'branch_scale'})/eg
		#my $prev = $_[0]; 
		#while ($_[0] =~ s/:(\d+\.\d+)/":" . ($1*$opts{'branch_scale'})/eg) {
		#	$prev eq $_[0] and last; 
		#	$prev = $_[0]; 
		#}
	}
	return $_[0]; 
}
