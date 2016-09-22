#!/usr/bin/perl -w
use strict; 
use warnings; 
use LogInforSunhh; 
use mathSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"bin_len:i", 
); 
my %default_para = ( 'bin_len' => 50 ); 

my $help_txt = <<HH; 
####################################################################################################
# perl $0 -bin_len 50 LD_grp1.LD_cnt > LD_grp1.LD_cnt_bin50bp 
#
# -bin_len       [$default_para{'bin_len'}]
####################################################################################################
### Input file format : 
# head -4 LD_grp1.LD_cnt
# Dist	pairNum	sum_R2	sum_Dp	Avg_r2	Avg_Dp
# 5	17960	6790.12	16622.716	0.378069042316258	0.925540979955457
# 6	16311	6227.476	15038.698	0.381796088529213	0.921997302433941
# 7	12913	4707.891	11803.864	0.36458537907535	0.914107023929374
####################################################################################################
HH

defined $opts{'help'} and &LogInforSunhh::usage($help_txt); 
-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

my %glob; 
&mathSunhh::_addHash( 'toH' => \%glob, 'fromH'=>\%opts, 'replaceExist' => 1 ); 
&mathSunhh::_addHash( 'toH' => \%glob, 'fromH'=>\%default_para, 'replaceExist' => 0 ); 

$glob{'_has_header'} = 0; 
my %cnt; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ( $ta[0] eq 'Dist' ) {
		$glob{'_has_header'} == 0 and do { print STDOUT "$_\n"; $glob{'_has_header'} = 1 }; 
		next; 
	}
	my $binN = int( ($ta[0]-1)/$glob{'bin_len'} ); 
	my $lvlN = int( $binN * $glob{'bin_len'} ) + 1 ; 
	$cnt{$lvlN}{'num'}    += $ta[1]; 
	$cnt{$lvlN}{'sum_R2'} += $ta[2]; 
	$cnt{$lvlN}{'sum_Dp'} += $ta[3]; 
}#

for my $lvlN (sort { $a <=> $b } keys %cnt) {
	print STDOUT join("\t", 
	  $lvlN, 
	  $cnt{$lvlN}{'num'}, 
	  $cnt{$lvlN}{'sum_R2'}, 
	  $cnt{$lvlN}{'sum_Dp'}, 
	  $cnt{$lvlN}{'sum_R2'} / $cnt{$lvlN}{'num'}, 
	  $cnt{$lvlN}{'sum_Dp'} / $cnt{$lvlN}{'num'}
	)."\n"; 
}

