#!/usr/bin/perl 
use strict; 
use warnings; 
use ReadInAlnSunhh; 
use LogInforSunhh; 
use Getopt::Long; 

my %opts; 
GetOptions(\%opts, 
	"help!", 
	"out:s", 
); 

-t and !@ARGV and &stopErr("[Err] fasdfaf\n"); 

my $oFh = \*STDOUT;
if (defined $opts{out}) {
	my $tfh;
	open $tfh,'>',"$opts{out}" or &stopErr("[Err] Failed to open $opts{out}\n$!\n");
	$oFh = $tfh;
}

my @FH;
!(-t) and push(@FH, \*STDIN);
for (@ARGV) {
	my $fh;
	open $fh, '<', "$_" or &stopErr("[Err] Failed to open $_\n$!\n");
	push(@FH, $fh);
}


my @all_blks;
my %links; 
&tsmsg("[Rec] Reading file.\n"); 
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
			# $specs[$i] ne '' and $sline{seqId} = "$specs[$i].$sline{seqId}";
			push(@cur_blk, [@sline{qw/seqId seqStart blkSize seqStrand seqLen seqSeq/}]);
			# print {$oFh} join(" ", "s", @sline{qw/seqId seqStart blkSize seqStrand seqLen seqSeq/})."\n";
		}
		$links{ $cur_blk[1][0] }{ $cur_blk[2][0] } = 1; 
		$links{ $cur_blk[2][0] }{ $cur_blk[1][0] } = 1; 
		push(@all_blks, [@cur_blk]);
		# print {$oFh} "\n";
	}
}
&tsmsg("[Rec] Files read ok.\n"); 

# Use %links to group scaffolds. 

&tsmsg("[Rec] Begin to group scaffolds\n"); 
my $grp_ct=0; 
my %used; 
for my $r1 (sort { $b->[1][4]<=>$a->[1][4] || $a->[1][3] cmp $b->[1][3] || $a->[1][1] <=> $b->[1][1]  } @all_blks) {
	defined $used{ $r1->[1][0] } and next; 
	my @grp_scf = &relatScaf( $r1->[1][0], \%links, {} ); 
	$grp_ct ++; 
	my ($refN, $qryN) = (0, 0); 
	for my $tmp_scf ( @grp_scf ) {
		defined $used{ $tmp_scf } and &tsmsg("[Err] Error because [$tmp_scf] has been used.\n"); 
		$tmp_scf =~ m/^Ref\./ and $refN ++; 
		$tmp_scf =~ m/^Qry\./ and $qryN ++; 
		$used{ $tmp_scf } = 1; 
	}
	print {$oFh} join("\t", "Group_$grp_ct", $refN, $qryN, join(" ", @grp_scf))."\n"; 
}
&tsmsg("[Rec] All done.\n"); 


sub relatScaf {
	my ( $base_scf, $linkR, $ignoreR ) = @_; 
	defined $ignoreR or $ignoreR = {}; 
	ref($ignoreR) eq 'HASH' or &stopErr("[Err] ignoreR is not HASH reference. [", ref($ignoreR), "]\n"); 
	my @rel_scf = ( $base_scf ); 
	$ignoreR->{ $base_scf } = 1; 
	if ( defined $linkR->{$base_scf} ) {
		my @rel_scf_2 = keys %{ $linkR->{$base_scf} }; 
		my @rel_scf_2_filt; 
		for my $scf2 ( @rel_scf_2 ) {
			defined $ignoreR->{ $scf2 } and next; 
			push(@rel_scf_2_filt, $scf2); 
			$ignoreR->{ $scf2 } = 1; 
		}
		for my $scf2 ( @rel_scf_2_filt ) {
			my @rel_scf_3 = &relatScaf( $scf2, $linkR, $ignoreR ); 
			push(@rel_scf, @rel_scf_3); 
		}
	} else {
		; 
	}

	return (@rel_scf); 
}


