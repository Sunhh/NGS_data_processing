#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"gm_id_list:s", 
); 

my (@mrk_list, $mrk_idx); 
$mrk_idx = -1; 
if (defined $opts{'gm_id_list'}) {
	@mrk_list = &fileSunhh::load_tabFile( $opts{'gm_id_list'}, 1 ); 
}

-t and !@ARGV and die "perl $0 in.snp.abh > in.snp.abh.rqtl.csv\n -gm_id_list in.snp.abh.mrk_id_list\n"; 


my %good = map { $_=>1 } qw/A B H -/; 

my $hl = <>; chomp($hl); 
$mrk_idx ++; 
my @ha = split(/\t/, $hl); 
my @out; 
$out[0] = ["id"]; 
$out[1] = [""]; 
$out[2] = [""]; 
for (my $i=0; $i<@ha-5; $i++) {
	$out[$i+3][0] = $ha[$i+5]; 
}
while (<>) {
	$. % 1000 == 1 and &tsmsg("[Msg] $. line\n"); 
	$mrk_idx ++; 
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[2] =~ s![^a-zA-Z0-9_]!_!g; 
	$ta[3] =~ s![^a-zA-Z0-9_]!_!g; 
	$ta[4] =~ s![^a-zA-Z0-9_]!_!g; 
	my $n = $ta[0]; 
	$n =~ s!^chr0*!!i; 
	$n eq '' and $n = 0; 
	my $p = 0; 
	my $snpID = "$ta[0]_$ta[3]_$ta[4]_$ta[2]"; 
	if ( defined $opts{'gm_id_list'} ) {
		$snpID = $mrk_list[$mrk_idx][0]; 
		$n = $mrk_list[$mrk_idx][1]; 
		$p = $mrk_list[$mrk_idx][2]; 
	}
	push(@{$out[0]}, $snpID); 
	push(@{$out[1]}, $n); 
	push(@{$out[2]}, $p); 
	for (my $i=5; $i<@ta; $i++) {
		( $ta[$i] eq 'u' or $ta[$i] eq 'U' ) and $ta[$i] = '-'; 
		$ta[$i] = uc($ta[$i]); 
		defined $good{$ta[$i]} or die "$ta[$i]\n"; 
		push(@{$out[$i-2]}, $ta[$i]); 
	}
}

for my $t1 (@out) {
	print join(",", @$t1)."\n"; 
}

