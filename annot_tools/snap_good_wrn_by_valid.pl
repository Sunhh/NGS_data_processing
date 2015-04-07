#!/usr/bin/perl
# Read in genome.fat_valid file, and find good models. 
# Read in wrn.dna and wrn.ann files, retrieve good models. 
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts,
	"valid:s", # genome.fat_valid
	"wrnDna:s", # wrn.dna
	"wrnAnn:s", # wrn.ann
	"outPref:s", # wrnGood
	"yaN:i", # 0123, similar to spaln -yaN 
	"help!", 
); 

$opts{'valid'} //= 'genome.fat_valid'; 
$opts{'wrnDna'} //= 'wrn.dna'; 
$opts{'wrnAnn'} //= 'wrn.ann'; 
$opts{'outPref'} //= 'wrnGood'; 
$opts{'yaN'} //= 2; 

sub usage {
print <<HH;
################################################################################
# perl $0 
#
# -help
#
# -valid         genome.fat_valid
# -wrnDna        wrn.dna
# -wrnAnn        wrn.ann
# -outPref       wrnGood . This will output two files: wrnGood.dna and wrnGood.ann
# -yaN           [2]
################################################################################
HH
	exit 1; 
}
$opts{'help'} and &usage(); 


open V,'<',"$opts{'valid'}" or &stopErr("[Err] -valid $!\n"); 
my %goodModel; 
while (<V>) {
	chomp; 
	if ( m/^([^:]+):\s+(\S+)\s+OK/ ) {
		# validated as good
		$goodModel{$2} = 1; 
	} elsif ( s!^([^:]+):\s+(\S+)\s+(\d+)\s+(\d+)\s+.*\swarnings\(\d+\):\s+!! ) {
		my $scfID = $1; 
		my $genID = $2; 
		my ($genS, $genE) = ($3, $4); 
		my @types = ($_ =~ m!(\S+)!g); 
		my $is_good = 1; 
		CHK_WRN: for my $tmp (@types) {
			if ( $tmp =~ m!^cds:incomplete$! ) {
				$is_good = 0; last CHK_WRN; 
			} elsif ( $tmp =~ m!^exon\-\d+:short\(\d+\)$! ) {
				; 
			} elsif ( $tmp =~ m!^split\-stop$! ) {
				$is_good = 0; last CHK_WRN; 
			} elsif ( $tmp =~ m!^intron\-\d+:(\w{2})\.\.(\w{2})$! ) {
				my $intron_seg1 = $1; 
				my $intron_seg2 = $2; 
				my $intron_seg12 = "${intron_seg1}..${intron_seg2}"; 
				if ( $opts{'yaN'} eq '0' and $intron_seg12 =~ m!^(GT\.\.AG|GC\.\.AG|AT\.\.AC)$! ) {
					; 
				} elsif ( $opts{'yaN'} eq '1' and $intron_seg12 =~ m!^(GT\.\.AG|GC\.\.AG|AT\.\.AC|AT\.\.A.)$! ) {
					; 
				} elsif ( $opts{'yaN'} eq '2' and $intron_seg12 =~ m!^(GT\.\.AG|GC\.\.AG|AT\.\.AC|AT\.\.A.|.T\.\.AG|G.\.\.AG|GT\.\..G|GT\.\.A.)$! ) {
					; 
				} elsif ( $opts{'yaN'} eq '3' ) {
					; 
				} else {
					$is_good = 0; 
					last CHK_WRN; 
				}
			} else {
				&stopErr("[Err] Unknown wrn_type: [$tmp]\n"); 
			}
		}#End CHK_WRN: 
		if ( $is_good == 1 ) {
			$goodModel{$genID} = 1; 
		}
	} elsif ( m!^E(init|xon|term)\t! ) {
		; 
	} else {
		&stopErr("[Err] Unknown line:$_\n"); 
	}
}
close V; 

open WA,'<',"$opts{'wrnAnn'}" or &stopErr("[Err] -wrnAnn $!\n"); 
my %infor_scfID; 
my ($scfID);
while (<WA>) {
	if (m!^>(\S+)!) {
		$scfID = $1; 
		$infor_scfID{$scfID}{'use'} = 0; 
		$infor_scfID{$scfID}{'lines'} .= $_; 
		$infor_scfID{$scfID}{'ln'} = $.; 
	}else{
		my @ta = split(/\t/, $_); 
		$ta[3] =~ s!\s+$!!; 
		if ( defined $goodModel{$ta[3]} and $goodModel{$ta[3]} == 1 ) {
			$infor_scfID{$scfID}{'use'} = 1; 
			$infor_scfID{$scfID}{'lines'} .= $_; 
		}
	}
}
close WA; 
open WGA,'>',"$opts{'outPref'}.ann" or &stopErr("[Err] -outPref $!\n"); 
for my $id ( sort { $infor_scfID{$a}{'ln'} <=> $infor_scfID{$b}{'ln'} } keys %infor_scfID ) {
	$infor_scfID{$id}{'use'} == 1 or next; 
	print WGA "$infor_scfID{$id}{'lines'}"; 
}
close WGA; 

open WD,'<',"$opts{'wrnDna'}" or &stopErr("[Err] -wrnDna $!\n"); 
open WGD,'>',"$opts{'outPref'}.dna" or &stopErr("[Err] -outPref $!\n"); 
my $is_out = 0; 
while (<WD>) {
	if (m!^\s*>(\S+)!) {
		if ( defined $infor_scfID{$1} and $infor_scfID{$1}{'use'} == 1 ) {
			$is_out = 1; 
			print WGD "$_"; 
		} else {
			$is_out = 0; 
		}
	} elsif ( $is_out == 1 ) {
		print WGD "$_"; 
	}
}
close WGD; 
close WD; 




