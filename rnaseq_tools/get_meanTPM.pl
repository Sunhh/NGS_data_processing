#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"list_sample2grp:s", 
); 

my %gg; 

&setup_Glob(); 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my @outLine = @ta[@{$gg{'outs'}}]; 
	for my $g0 (@{$gg{'grps'}}) {
		my ($sum, $cnt)=(0,0); 
		for my $i (@{$g0->[1]}) {
			$sum += $ta[$i]; 
			$cnt ++; 
		}
		push(@outLine, $sum/$cnt); 
	}
	print STDOUT join("\t", @outLine)."\n"; 
}

sub setup_Glob {
	$gg{'h_txt'} = <<"H1"; 
################################################################################
# perl $0   rdCnt.Both.TPM   > rdCnt.Both.TPMmean
#
# -help 
#
# -list_sample2grp        [filename] Only samples listed in this file will be grouped if given. 
#                           Foramt: sampleID1 \\t grpID1
#                                   sampleID2 \\t grpID1
#                                   sampleID3 \\t grpID2
#                                   sampleID1 \\t grpID2
#                                   .... 
################################################################################
H1
	
	-t and !@ARGV and &LogInforSunhh::usage($gg{'h_txt'}); 
	$opts{'help'} and &LogInforSunhh::usage($gg{'h_txt'}); 
	
	$gg{'hl_txt'} = <>; chomp($gg{'hl_txt'}); $gg{'hl_arr'} = [ &splitL("\t", $gg{'hl_txt'}) ]; 
	if (defined $opts{'list_sample2grp'}) {
		my $fh = &openFH($opts{'list_sample2grp'}, '<'); 
		while (<$fh>) {
			m!^\s*($|#)! and next; 
			chomp; 
			my @ta = &splitL("\t", $_); 
			push(@{$gg{'lis_s2g'}{$ta[0]}}, $ta[1]); 
		}
		close($fh); 
	}
	$gg{'grps'} = []; 
	$gg{'outs'} = []; 
	for (my $i=0; $i<@{$gg{'hl_arr'}}; $i++) {
		my @arr_grpID; 
		if ( defined $opts{'list_sample2grp'} ) {
			if ( defined $gg{'lis_s2g'}{$gg{'hl_arr'}[$i]} ) {
				@arr_grpID = @{ $gg{'lis_s2g'}{$gg{'hl_arr'}[$i]} }; 
			} elsif ( $i == 0 or $gg{'hl_arr'}[$i] =~ m!^(eleID|geneID|mrnaID)$!i ) {
				&tsmsg("[Wrn] Keep column of [$gg{'hl_arr'}[$i]]\n"); 
				push(@{$gg{'outs'}}, $i); 
				next; 
			} else {
				next; 
			}
		} else {
			$gg{'hl_arr'}[$i] =~ m!^(\S+)_Rep(\d+)$! or do { &tsmsg("[Wrn] Skip bad ID [$gg{'hl_arr'}[$i]]\n"); push(@{$gg{'outs'}}, $i); next;  }; 
			@arr_grpID = ($1); 
		}
		for my $grpID (@arr_grpID) {
			if ( @{$gg{'grps'}} == 0 ) {
				push(@{$gg{'grps'}}, [$grpID, [$i], { $gg{'hl_arr'}[$i]=>1 }]); 
				next; 
			}
			my $is_foundGrp = 0; 
			for my $g0 (@{$gg{'grps'}}) {
				if ($g0->[0] eq $grpID) {
					$is_foundGrp = 1; 
					if ( defined $g0->[2]{$gg{'hl_arr'}[$i]} ) {
						&tsmsg("[Wrn] Skip repeated sample ID [$gg{'hl_arr'}[$i]]\n"); 
					} else {
						push(@{$g0->[1]}, $i); 
						$g0->[2]{$gg{'hl_arr'}[$i]} = 1; 
					}
				}
			}
			$is_foundGrp == 1 and last; 
			if ($is_foundGrp == 0) {
				push(@{$gg{'grps'}}, [$grpID, [$i], { $gg{'hl_arr'}[$i]=>1 }]); 
			}
		}# End for my $grpID 
	}
	
	@{$gg{'grps'}} > 0 or &stopErr("[Err] No good sample ID found!\n"); 
	my @outHeader = @{$gg{'hl_arr'}}[@{$gg{'outs'}}]; 
	for my $g0 (@{$gg{'grps'}}) {
		my $grpID = "$g0->[0]_x" . scalar(@{$g0->[1]}); 
		push(@outHeader, $grpID); 
	}
	print STDOUT join("\t", @outHeader)."\n"; 
	# 
}# setup_Glob() 

