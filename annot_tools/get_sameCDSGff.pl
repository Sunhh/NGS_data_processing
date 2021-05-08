#!/usr/bin/perl
### Only CDS blocks are considered. 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"srcGff:s", 
	"idxGff:s", 
	"exact!", 
	"keepCDSorder!", 
); 

my $help_txt = <<HH; 
######################################################################
# perl $0 -srcGff maker_r1.gff3 -idxGff maker_r2.gff3 > maker_r1_same2r2.gff3
######################################################################
###  Only feature IDs "mRNA" and "CDS" are considered. 
#
# -srcGff     [Filename] 
# -idxGff     [Filename] as index
#
# -exact      [Boolean] Require the first and last CDS block has the same boundries too. 
# 
# -keepCDSorder [Boolean] Keep the order of CDS blocks as input. 
HH

( defined $opts{'srcGff'} and defined $opts{'idxGff'} ) or &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 


my %trans_1 = %{ &load_gff3($opts{'srcGff'}) }; 
my %trans_2 = %{ &load_gff3($opts{'idxGff'}) }; 
my %cIntron_1 = %{ &hash_CDSintron(\%trans_1) }; 
my %cIntron_2 = %{ &hash_CDSintron(\%trans_2) }; 
my %cList_1; 
my %cList_2; 
%cList_1 = %{ &hash_CDSlist(\%trans_1) }; 
%cList_2 = %{ &hash_CDSlist(\%trans_2) }; 

for my $mID_1 ( sort keys %trans_1 ) {
	defined $trans_1{$mID_1}{'cds'} or next; 
	defined $trans_1{$mID_1}{'CDSintron'} or &stopErr("[Err] Failed to find CDSintron for [$mID_1]\n"); 
	my $cIntron_test1 = $trans_1{$mID_1}{'CDSintron'}; 
	defined $cIntron_2{$cIntron_test1} or next; 
	if ( $opts{'exact'} ) {
		my $cList_test1 = $trans_1{$mID_1}{'CDSlist'}; 
		defined $cList_2{$cList_test1} or next; 
	}
	print STDOUT $trans_1{$mID_1}{'txt'}; 
}

sub load_gff3 {
	my ($fn) = @_; 
	&tsmsg("[Msg] Loading gff3 [$fn]\n"); 
	my $fh = &openFH($fn) or &stopErr("[Err] Failed to load [$fn]\n"); 
	my %transInfo; 
	my %tmp_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>1e5 ); 
	while (<$fh>) {
		&fileSunhh::log_section( $., \%tmp_cnt ) and &tsmsg("[Msg]   Reading $. line.\n"); ; 
		m!^##FASTA\s*$! and last; 
		m!^\s*(#|$)! and next; 
		chomp; 
		my @ta = split(/\t/, $_); 
		my $fname = lc($ta[2]); 
		if ( $fname eq 'mrna' ) {
			$ta[8] =~ m!(?:^|\s|;)ID=([^\s;]+)\s*(?:;|$)! or die "bad ta[8] 1 [$ta[8]]\n"; 
			my $mID = $1; 
			$mID =~ m!,! and &stopErr("[Err] No ',' allowed in ID; [$mID]\n"); 
			if (defined $transInfo{$mID} and defined $transInfo{$mID}{'SE'}) {
				&stopErr("[Err] Repeat ID for mID [$mID]\n"); 
			}
			$transInfo{$mID}{'str'} = $ta[6]; 
			$transInfo{$mID}{'SE'} = [@ta[0,3,4]]; 
			$transInfo{$mID}{'txt'} .= "$_\n"; 
		} elsif ( $fname eq 'cds' ) {
			$ta[8] =~ m!(?:^|\s|;)ID=([^\s;]+)\s*(?:;|$)! or die "bad ta[8] 2 [$ta[8]]\n"; 
			my $cID = $1; 
			$cID =~ m!,! and &stopErr("[Err] No ',' allowed in ID; [$cID]\n");
			$ta[8] =~ m!(?:^|\s|;)Parent=([^\s;]+)\s*(?:;|$)! or die "bad ta[8] 3 [$ta[8]]\n"; 
			my $pIDs = $1; 
			for my $pID (split(/,/,$pIDs)) {
				$pID =~ s!^\s+|\s+$!!g; 
				$transInfo{$pID}{'str'} //= $ta[6]; 
				if (defined $transInfo{$pID}{'cds'}) {
					$transInfo{$pID}{'cds'}[-1][0] ne $ta[0] and !($opts{'keepCDSorder'}) and &stopErr("[Err] Please use -keepCDSorder for gene case $pID\n"); 
				}
				push(@{$transInfo{$pID}{'cds'}}, [@ta[0,3,4]]); 
				$transInfo{$pID}{'txt'} .= "$_\n"; 
			}
		} else {
			next; 
		}
	}
	close($fh); 
	if ( !($opts{'keepCDSorder'}) ) {
		for my $mID (keys %transInfo) {
			defined $transInfo{$mID}{'cds'} or next; 
			( defined $transInfo{$mID}{'str'} and $transInfo{$mID}{'str'} =~ m!^(\+|\-)$! ) or die "str=[$transInfo{$mID}{'str'}] $mID\n"; 
			if ( $transInfo{$mID}{'str'} eq '+' ) {
				@{$transInfo{$mID}{'cds'}} = sort { $a->[1] <=> $b->[1] || $a->[1] <=> $b->[2] } @{$transInfo{$mID}{'cds'}}; 
			} elsif ( $transInfo{$mID}{'str'} eq '-' ) {
				@{$transInfo{$mID}{'cds'}} = sort { $b->[1] <=> $a->[1] || $b->[1] <=> $a->[2] } @{$transInfo{$mID}{'cds'}}; 
			}
		}
	}
	for my $mID (keys %transInfo) {
		defined $transInfo{$mID}{'cds'} or next;
		my @ii; 
		for (my $i=1; $i<@{$transInfo{$mID}{'cds'}}; $i++) {
			defined $transInfo{$mID}{'cds'} or next; 
			my $is = $transInfo{$mID}{'cds'}[$i-1][2]+1; $is = "$transInfo{$mID}{'cds'}[$i-1][0]_+_$is"; 
			my $ie = $transInfo{$mID}{'cds'}[$i  ][1]-1; $ie = "$transInfo{$mID}{'cds'}[$i  ][0]_+_$ie"; 
			if ($transInfo{$mID}{'str'} eq '-') {
				$is = $transInfo{$mID}{'cds'}[$i-1][1]-1; $is = "$transInfo{$mID}{'cds'}[$i-1][0]_-_$is";
				$ie = $transInfo{$mID}{'cds'}[$i  ][2]+1; $ie = "$transInfo{$mID}{'cds'}[$i  ][0]_-_$ie";
			}
			push(@ii, "$is-$ie"); 
		}
		if (@ii == 0) {
			if ($transInfo{$mID}{'str'} eq '+') {
				$transInfo{$mID}{'CDSintron'} = "$transInfo{$mID}{'cds'}[0][0]_+0_$transInfo{$mID}{'cds'}[0][1]-$transInfo{$mID}{'cds'}[0][0]_+0_$transInfo{$mID}{'cds'}[0][2]"; 
			} else {
				$transInfo{$mID}{'CDSintron'} = "$transInfo{$mID}{'cds'}[0][0]_-0_$transInfo{$mID}{'cds'}[0][2]-$transInfo{$mID}{'cds'}[0][0]_-0_$transInfo{$mID}{'cds'}[0][1]"; 
			}
		} else {
			$transInfo{$mID}{'CDSintron'} = join(";", @ii) ; 
		}
		$transInfo{$mID}{'CDSlist'}   = join(";", map { "$_->[0]\t$transInfo{$mID}{'str'}\t$_->[1]\t$_->[2]" } @{$transInfo{$mID}{'cds'}}); 
	}
	return(\%transInfo); 
}# load_gff3() 

sub hash_CDSintron {
	my ($transH) = @_; 
	my %h; 
	for my $mID (keys %$transH) {
		( defined $transH->{$mID}{'CDSintron'} and $transH->{$mID}{'CDSintron'} ne '' ) or next; 
		$h{ $transH->{$mID}{'CDSintron'} } ++; 
	}
	return(\%h); 
}# hash_CDSintron() 
sub hash_CDSlist {
	my ($transH) = @_; 
	my %h; 
	for my $mID (keys %$transH) {
		( defined $transH->{$mID}{'CDSlist'} and $transH->{$mID}{'CDSlist'} ne '' ) or next; 
		$h{ $transH->{$mID}{'CDSlist'} } ++; 
	}
	return(\%h); 
}# hash_CDSlist() 


