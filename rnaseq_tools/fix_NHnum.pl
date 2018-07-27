#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use SeqAlnSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"inBam:s", # 
	"outBam:s", 
	"nSorted!", 
	"para_samTsrt:s", 
	"exe_samtools:s", 
	"help!", 
); 

$opts{'exe_samtools'} //= 'samtools'; 
$opts{'para_samTsrt'} //= ' -@ 4 -m 5G '; 

my $help_txt = <<HH; 
################################################################################
# perl $0   -inBam hisat2_filtered.bam -outBam hisat2_fixNH.bam 
#
# -help 
#
# -exe_samtools   ['$opts{'exe_samtools'}']
# -para_samTsrt   ['$opts{'para_samTsrt'}']
# -nSorted        [Boolean] Already sorted by name if given. 
################################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
defined $opts{'inBam'}  or &LogInforSunhh::usage($help_txt); 
defined $opts{'outBam'} or &LogInforSunhh::usage($help_txt); 


my %flag_RS = %{ &SeqAlnSunhh::mk_flag( 'keep'=>'0=0,6=0,7=0' ) }; 
my %flag_R1 = %{ &SeqAlnSunhh::mk_flag( 'keep'=>'0=1,6=1,7=0' ) }; 
my %flag_R2 = %{ &SeqAlnSunhh::mk_flag( 'keep'=>'0=1,6=0,7=1' ) }; 
my %flag2R; 
for (keys %flag_RS) { $flag2R{$_} = 'RS'; }
for (keys %flag_R1) { $flag2R{$_} = 'R1'; }
for (keys %flag_R2) { $flag2R{$_} = 'R2'; }

my $tmpDir = &fileSunhh::new_tmp_dir('create'=>1); 

my $baseFn = &fileSunhh::_basename( $opts{'inBam'} ); 
if ( $opts{'nSorted'} ) {
	open F,'-|', "$opts{'exe_samtools'} view -h $opts{'inBam'}" or die "$!\n"; 
} else {
	&exeCmd_1cmd("$opts{'exe_samtools'} sort -n $opts{'para_samTsrt'} -o $tmpDir/nSrt_$baseFn $opts{'inBam'}"); 
	open F,'-|', "$opts{'exe_samtools'} view -h $tmpDir/nSrt_$baseFn" or die "$!\n"; 
}
open O,'|-', "$opts{'exe_samtools'} view -o $opts{'outBam'} -" or die "$!\n"; 
my @sameRd=('', [], 0); 
while (<F>) {
	if (m!^\@!) {
		print O; 
		next; 
	}
	chomp; 
	my @ta=split(/\t/, $_); 
	my $r = $flag2R{$ta[1]}; 
	( defined $r and $r ne '' ) or die "Bad flag [$ta[1]]\n"; 
	my $rdKey = "$ta[0]\t$r"; 
	for (my $i=11; $i<@ta; $i++ ) {
		$ta[$i] =~ m!^NH:i:\d+$! or next; 
		splice(@ta, $i, 1); 
		last; 
	}
	if ($sameRd[0] eq $rdKey) {
		push(@{$sameRd[1]}, join("\t", @ta)); 
		$ta[2] ne '*' and $sameRd[2] ++; 
	} else {
		if ( $sameRd[2] > 0 ) {
			for my $tline (@{$sameRd[1]}) {
				print O "$tline\tNH:i:$sameRd[2]\n"; 
			}
		} else {
			for my $tline (@{$sameRd[1]}) {
				print O "$tline\n"; 
			}
		}
		@sameRd = ( $rdKey, [ join("\t", @ta) ], 0); 
		$ta[2] ne '*' and $sameRd[2] ++; 
	}
}
if ( $sameRd[2] > 0 ) {
	for my $tline (@{$sameRd[1]}) {
		print O "$tline\tNH:i:$sameRd[2]\n"; 
	}
} else {
	for my $tline (@{$sameRd[1]}) {
		print O "$tline\n"; 
	}
}
@sameRd = ( '', [ ], 0); 
close O; 
close F; 

&fileSunhh::_rmtree($tmpDir); 

&tsmsg("[Rec] Finished $opts{'inBam'}\n"); 

