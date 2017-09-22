#!/usr/bin/perl
use strict; 
use warnings; 
use ReadInAlnSunhh; 
use LogInforSunhh; 
use Getopt::Long;

my %opts; 
GetOptions(\%opts, 
	"help!", 
	"specs:s", 
	"out:s", 
	"minLen1:i", 
	"minLen2:i", 
); 

-t and !@ARGV and &usage(); 
$opts{help} and &usage(); 

sub usage {
print STDOUT <<HELP; 
perl $0 in.maf
-help 
-specs      [S01,S02] Name of species
-out        [\\*STDOUT] Output file name. 
-minLen1    [1] 
-minLen2    [1]
HELP
	exit 1; 
}

my @specs = qw(S01 S02); 
defined $opts{specs} and @specs = split(',', $opts{specs}); 
defined $specs[0] or $specs[0] = ''; 
defined $specs[1] or $specs[1] = ''; 
my $oFh = \*STDOUT; 
if (defined $opts{out}) {
	my $tfh; 
	open $tfh,'>',"$opts{out}" or &stopErr("[Err] Failed to open $opts{out}\n$!\n"); 
	$oFh = $tfh; 
}
defined $opts{minLen1} or $opts{minLen1} = 1; 
defined $opts{minLen2} or $opts{minLen2} = 1; 

my @FH; 
!(-t) and push(@FH, \*STDIN); 
for (@ARGV) {
	my $fh; 
	open $fh, '<', "$_" or &stopErr("[Err] Failed to open $_\n$!\n"); 
	push(@FH, $fh); 
}

print {$oFh} "##maf version=1\n"; 
my @all_blks; 
for my $fh (@FH) {
	while ( my %rec1 = %{readMAF($fh)} ) {
		my @cur_blk; 
		chomp( $rec1{a}[0] ); 
		push(@cur_blk, [ $rec1{a}[0] ]); 
		# print {$oFh} "$rec1{a}[0]\n"; 
		@{$rec1{o}} >= 2 or next; 
		for (my $i=0; $i<2; $i++) {
			$rec1{o}[$i] =~ m/^s\s/ or &stopErr("Wrong line[$i]: $rec1{o}[$i]\n"); 
			my %sline = %{ splitMafSline($rec1{o}[$i], 1)}; 
			$specs[$i] ne '' and $sline{seqId} = "$specs[$i].$sline{seqId}"; 
			push(@cur_blk, [@sline{qw/seqId seqStart blkSize seqStrand seqLen seqSeq/}]); 
			# print {$oFh} join(" ", "s", @sline{qw/seqId seqStart blkSize seqStrand seqLen seqSeq/})."\n"; 
		}

		$cur_blk[1][4] >= $opts{minLen1} or next; 
		$cur_blk[2][4] >= $opts{minLen2} or next; 

		push(@all_blks, [@cur_blk]); 
		# print {$oFh} "\n"; 
	}
}

for my $r1 (sort { $b->[1][4]<=>$a->[1][4] || $a->[1][3] cmp $b->[1][3] || $a->[1][1] <=> $b->[1][1]  } @all_blks) {
	print {$oFh} $r1->[0][0] . "\n"; 
	for (my $i=1; $i<@$r1; $i++) {
		print {$oFh} join(" ", "s", @{$r1->[$i]})."\n"; 
	}
	print {$oFh} "\n"; 
}


