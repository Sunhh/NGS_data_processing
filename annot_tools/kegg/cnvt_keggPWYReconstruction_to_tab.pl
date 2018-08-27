#!/usr/bin/perl
# 20180824 : The input file, 'keggPWYByKO_html2text.txt', comes from "KEGG Mapper – Reconstruct Pathway" tool : https://www.kegg.jp/kegg/tool/map_pathway.html . 
#
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"print_level_0!", 
	"load_level_0:s", 
	"load_level_skip:s", 
); 

my %gg; 
&set_Glob(); 
sub set_Glob {
	$gg{'level_N'} = -1; 
	$gg{'level_skip'}{'Pathway Reconstruction Result'} = 1; 
	$gg{'level_skip'}{'Hide all objects'}              = 2; 
	$gg{'level_0'}{'Metabolism'}                           = 1; 
	$gg{'level_0'}{'Genetic Information Processing'}       = 2; 
	$gg{'level_0'}{'Environmental Information Processing'} = 3; 
	$gg{'level_0'}{'Cellular Processes'}                   = 4; 
	$gg{'level_0'}{'Organismal Systems'}                   = 5; 
	$gg{'level_0'}{'Human Diseases'}                       = 6; 
	$gg{'htxt'} = <<"H1"; 
################################################################################
# perl $0 keggPWYByKO_html2text.txt > keggPWYByKO.tab
#
# -help
#
# -print_level_0        [Boolean] If given, print names of level_0 and exit. 
# -load_level_0         [filename] If given, default names of level_0 will be replaced. 
#                        Format : Metabolism             \\t 1
#                                 Cellular Processes     \\t 4
#                                 ....
# -load_level_skip      [filename] similar to -load_level_0
################################################################################
H1
	if ($opts{'print_level_0'}) {
		for my $id (sort keys %{$gg{'level_0'}}) {
			print STDOUT "$id\t$gg{'level_0'}{$id}\n"; 
		}
		exit(); 
	}
	!@ARGV and -t and &LogInforSunhh::usage($gg{'htxt'}); 
	$opts{'help'} and &LogInforSunhh::usage($gg{'htxt'}); 
	if (defined $opts{'load_level_0'}) {
		$gg{'level_0'} = {}; 
		my $n = 0; 
		for my $lr (&fileSunhh::load_tabFile($opts{'load_level_0'}, 1)) {
			$n ++; 
			$lr->[1] //= $n; 
			$gg{'level_0'}{ $lr->[0]     } = $lr->[1]; 
			$gg{'level_0'}{ lc($lr->[0]) } = $lr->[1]; 
		}
	}
	if (defined $opts{'load_level_skip'}) {
		$gg{'level_skip'} = {}; 
		my $n = 0; 
		for my $lr (&fileSunhh::load_tabFile($opts{'load_level_skip'}, 1)) {
			$n ++; 
			$lr->[1] //= $n; 
			$gg{'level_skip'}{ $lr->[0]     } = $lr->[1]; 
			$gg{'level_skip'}{ lc($lr->[0]) } = $lr->[1]; 
		}
	}
	for my $id (keys %{$gg{'level_0'}}) { $gg{'level_0'}{ lc($id) } = $gg{'level_0'}{ $id }; }
	for my $id (keys %{$gg{'level_skip'}}) { $gg{'level_skip'}{ lc($id) } = $gg{'level_skip'}{ $id }; }
}# set_Glob() 


# print STDOUT join("\t", qw/PWY_Group1 PWY_Group2 Map_ID Map_description KO_Num_InMap KO_ID Org_Idx Genes/)."\n"; 
while (<>) {
	chomp; 
	s/[\r\n]+$//; 
	if ( defined $gg{'level_skip'}{ lc($_) } ) {
		# Example : 'Pathway Reconstruction Result'
		$gg{'level_N'} = -1; 
		$gg{'levels'} = []; 
	} elsif (m!^\s*(\d+)\s+(\S.*)\s+\((\d+)\)$! and $gg{'level_N'} == 1) {
		# Example : '01100 Metabolic pathways (831)'
		my ($mapID, $mapDesc, $koNumInMap) = ($1, $2, $3); 
		@{$gg{'levels'}} = (@{$gg{'levels'}}[0,1], "$mapID\t$mapDesc\t$koNumInMap"); 
		$gg{'level_N'} = 2; 
	} elsif (m!^\s*(\d+)\s+(\S.*)\s+\((\d+)\)$! and $gg{'level_N'} >= 2) {
		my ($mapID, $mapDesc, $koNumInMap) = ($1, $2, $3); 
		@{$gg{'levels'}} = (@{$gg{'levels'}}[0,1], "$mapID\t$mapDesc\t$koNumInMap"); 
		$gg{'level_N'} = 2; 
	} elsif (m!^\s{2}(K\d+)\s*$! and $gg{'level_N'} == 2) {
		# Example : '  K00001'
		@{$gg{'levels'}} = (@{$gg{'levels'}}[0..2], $1); 
		$gg{'level_N'} = 3; 
		$gg{'orgNum'} = 0; 
	} elsif (m!^\s{2}(K\d+)\s*$! and $gg{'level_N'} >= 3) {
		@{$gg{'levels'}} = (@{$gg{'levels'}}[0..2], $1); 
		$gg{'level_N'} = 3; 
		$gg{'orgNum'} = 0; 
	} elsif (defined $gg{'level_0'}{$_}) {
		$gg{'levels'} = [$_]; 
		$gg{'level_N'} = 0; 
	} elsif (m!^(\S.*)$! and $gg{'level_N'} == 0) {
		# Example : 'Global and overview maps'
		@{$gg{'levels'}} = ($gg{'levels'}[0], $1); 
		$gg{'level_N'} = 1; 
	} elsif (m!^(\S.*)$! and $gg{'level_N'} >= 1) {
		@{$gg{'levels'}} = ($gg{'levels'}[0], $1); 
		$gg{'level_N'} = 1; 
	} elsif (m!^\s{3}(\S.*|)$! and $gg{'level_N'} == 3) {
		$gg{'orgNum'} ++; 
		$gg{'orgNumMax'} //= $gg{'orgNum'}; 
		$gg{'orgNumMax'} < $gg{'orgNum'} and $gg{'orgNumMax'} = $gg{'orgNum'}; 
		# print STDOUT join("\t", @{$gg{'levels'}}[0 .. 3], $gg{'orgNum'}, $1)."\n"; 
		push(@{$gg{'outSource'}}, [ @{$gg{'levels'}}[0 .. 3], $gg{'orgNum'}, "$1" ]); 
	} elsif (m!^\s*$!) {
		if ( $gg{'level_N'} == -1 ) {
			next; 
		} else {
			die "1:$.:$_\n"; 
		}
	} else {
		die "2:$_\nlevel_N=$gg{'level_N'}\n"; 
	}
}

print STDOUT join("\t", qw/PWY_Group1 PWY_Group2 Map_ID Map_description KO_Num_InMap KO_ID/, map { "Org_$_" } (1 .. $gg{'orgNumMax'}))."\n"; 
for (my $i=0; $i<@{$gg{'outSource'}}; $i+= $gg{'orgNumMax'}) {
	my @genes; 
	for (my $j=0; $j<$gg{'orgNumMax'}; $j++) {
		my $k = $i+$j; 
		$gg{'outSource'}[$k][4] == $j+1 or &stopErr("[Err] Unequal organism idx for line [$i+$j] org[$gg{'outSource'}[$k][4] : $j+1]" . join("\t", @{$gg{'outSource'}[$k]})."\n"); 
		push(@genes, $gg{'outSource'}[$k][5]); 
	}
	print STDOUT join("\t", @{$gg{'outSource'}[$i]}[0..3], @genes)."\n"; 
}


