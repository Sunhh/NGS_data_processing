#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use mathSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"processed!", 
	"tagPresent!", 
	"cntN_step:i", 
	"help!", 
); 
$opts{'cntN_step'} //= 100; 

my %tmp_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>$opts{'cntN_step'} ); 

my $help_txt = <<HH; 
######################################################################
### This is used to group SNP sites which have same distribution among input acessions. 
# gzip -cd WmReseq_sample473_rmClose_rmSame_TreeID.srtID.snp.gz | perl $0 > grped_lines\n
# 
# -processed        [Boolean] If given, the input file is treated as result of $0; 
# -tagPresent       [Boolean] If given, the output file will have '+' before representative lines; 
# -cntN_step        [Number] Default [$opts{'cntN_step'}]; 
#
# version 3.0
######################################################################
HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

# A-1 T-2 G-3 C-4 Ins-5 Del-6 N-0
# AA-11 TT-22 ...

my %c2v = qw(
 N   0
 A   1
 T   2
 C   3
 G   4
 Ins 5
 Del 6
 *   6
); 
my $n_geno = $c2v{'N'} . $c2v{'N'}; 
my $n_add  = 'AA'; 

my @g1_snp_line; 
# my @cur_fmt_line; 

my $h = <>; 
print STDOUT $h; 
while (<>) {
	if ( &fileSunhh::log_section($., \%tmp_cnt) ) {
		my $grpN = scalar(@g1_snp_line); 
		my $lN = $.-2; 
		&tsmsg("[Msg] Processed $lN line in $grpN groups.\n"); 
	}
	chomp; 
	my @ta = split(/\t/, $_); 
	my @ta1 = splice(@ta, 0, 2); 
	if ( $opts{'processed'} ) {
		my $tag_present = ''; 
		$ta1[0] =~ s!^(\+)!! and $tag_present = $1; 
		if ( !defined $ta[0] or $ta[0] eq '' ) {
			@g1_snp_line == 0 and next; 
			push(@{$g1_snp_line[-1][0]}, [@ta1] ); 
			$tag_present ne '' and $g1_snp_line[-1][4] = [@ta1]; 
			next; 
		}
		my %cur_snp; 
		my $cur_gIdx = -1; 
		for my $b (@ta) {
			$cur_gIdx ++; 
			push(@{$cur_snp{'grp_indv'}{$b}}, $cur_gIdx); 
		}
		@{$cur_snp{'key'}} = grep { $_ ne $n_geno } sort keys %{$cur_snp{'grp_indv'}}; 
		$cur_snp{'loc'} = [[@ta1]]; 
		$cur_snp{'line'} = [@ta]; 
		push(@g1_snp_line, 
			[
				[ @{$cur_snp{'loc'}} ], 
				[ @{$cur_snp{'line'}} ], 
				{ %{$cur_snp{'grp_indv'}} }, 
				[ @{$cur_snp{'key'}} ],
				[ @{$cur_snp{'loc'}[0]} ]
			]
		); 
		next; 
	}# processed; 
	# Format SNP chars to numbers, and group indv; 
	my %cur_snp; 
	my $cur_gIdx = -1; 
	for my $b (@ta) {
		$cur_gIdx ++; 
		$b = uc($b); 
		if ($b =~ m!^[ATGCN*]$!) {
			$b = $c2v{$b} . $c2v{$b}; 
		} elsif ($b =~ m!^([ATGCN*])([ATGCN*])$!) {
			$b = $c2v{$1} . $c2v{$2}; 
		} elsif ($b =~ m!^[ATGCN*]*\+[ATGCN]+$!) {
			$b = $c2v{'Ins'} . $c2v{'Ins'}; 
#		} elsif ($b =~ m!^([ATGCN*])[ATGCN*]\+[ATGCN]+$!) {
#			$b = $c2v{'Ins'} . $c2v{'Ins'}; 
#			# $b = $c2v{$1} . $c2v{'Ins'}; 
		} elsif ($b =~ m!^[ATGCN*]{3,}$!) {
			$b = $c2v{'N'} . $c2v{'N'}; 
		} else {
			die "|$b|\n$_\n"; 
		}
		push(@{$cur_snp{'grp_indv'}{$b}}, $cur_gIdx); 
	}
	@{$cur_snp{'key'}} = grep { $_ ne $n_geno and $_ ne $n_add } sort keys %{$cur_snp{'grp_indv'}}; 
	$cur_snp{'loc'} = [[$ta1[0], $ta1[1]]]; 
	$cur_snp{'line'} = [@ta]; 
	$cur_snp{'present'} = [@ta1]; 

	# Compare to stored groups; 
	my $inGrp = 0; 
	CHK_GRP1: 
	for my $g1 (@g1_snp_line) {
		$inGrp = &update_grp( $g1, \%cur_snp ); 
		$inGrp == 1 and last CHK_GRP1; 
	}# for my $g1 () 

	$inGrp == 1 and next; 
	push(@g1_snp_line, 
		[
			[ @{$cur_snp{'loc'}} ], 
			[ @{$cur_snp{'line'}} ], 
			{ %{$cur_snp{'grp_indv'}} }, 
			[ @{$cur_snp{'key'}} ],
			[ @{$cur_snp{'loc'}[0]} ]
		]
	); 
}

my $prev_grpN = scalar(@g1_snp_line);
my $curr_grpN = 0; 
my $round = 0; 
while ( $prev_grpN > $curr_grpN ) {
	# repeatly compress @g1_snp_line; 
	$prev_grpN = scalar(@g1_snp_line); 
	$round ++; 
	&tsmsg("[Msg] Running round $round [$prev_grpN]\n"); 
	my @grp2; 
	for my $g1 (@g1_snp_line) {
		my $inGrp = 0; 
		my %cur_snp; 
		$cur_snp{'loc'} = $g1->[0]; 
		$cur_snp{'line'} = $g1->[1]; 
		$cur_snp{'grp_indv'} = $g1->[2]; 
		$cur_snp{'key'} = $g1->[3]; 
		$cur_snp{'present'} = $g1->[4]; 
		CHK_GRP2:
		for my $g2 (@grp2) {
			$inGrp = &update_grp( $g2, \%cur_snp ); 
			$inGrp == 1 and last CHK_GRP2; 
		}
		$inGrp == 1 and next; 
		push(@grp2, 
			[
				[ @{$cur_snp{'loc'}} ], 
				[ @{$cur_snp{'line'}} ], 
				{ %{$cur_snp{'grp_indv'}} }, 
				[ @{$cur_snp{'key'}} ], 
				[ @{$cur_snp{'present'}} ]
			]
		); 
	}# for my $g1 
	@g1_snp_line = (); 
	@g1_snp_line = @grp2; 
	$curr_grpN = scalar(@g1_snp_line); 
	&tsmsg("[Msg] After round $round [$curr_grpN]\n"); 
}# End while

# Sort loci in @grp_snp_line
for my $g1 (@g1_snp_line) {
	@{$g1->[0]} = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } @{$g1->[0]}; 
}
@g1_snp_line =  sort { $a->[0][0][0] cmp $b->[0][0][0] || $a->[0][0][1] <=> $b->[0][0][1] } @g1_snp_line; 

for my $g1 (@g1_snp_line) {
	my $n=0; 
	for my $p1 (@{$g1->[0]}) {
		$n ++; 
		my $tag = ''; 
		if ($opts{'tagPresent'} and $g1->[4][0] eq $p1->[0] and $g1->[4][1] == $p1->[1]) {
			$tag = '+'; 
		}
		if ($n == 1) {
			print STDOUT join("\t", "$tag$p1->[0]", $p1->[1], @{$g1->[1]})."\n"; 
			next; 
		}
		print STDOUT join("\t", "$tag$p1->[0]", $p1->[1])."\n"; 
	}
}


###################################################################################
###################################################################################
###################################################################################
sub compare_grpLen {
	my ($rg1, $rg2) = @_; 
	my $l1 = scalar(@$rg1); 
	my $l2 = scalar(@$rg2); 
	if ($l1 == $l2) {
		return('same'); 
	} elsif ($l1 > $l2) {
		return('more'); 
	} elsif ($l1 < $l2) {
		return('less'); 
	} else {
		die "Bad :6:\n"; 
		return(); 
	}
}# compare_grpLen()
sub compare_grpContent {
	# Return (diff|same|more|less); # rg1 is XXXX than rg2; 
	my ($rg1, $rg2) = @_; 
	if (@$rg1 == @$rg2) {
		my %h1 = map { $_=>1 } @$rg1; 
		for (@$rg2) {
			defined $h1{$_} or return('diff'); 
		}
		return('same'); 
	} elsif (@$rg1 > @$rg2) {
		my %h1 = map { $_=>1 } @$rg1; 
		for (@$rg2) {
			defined $h1{$_} or return('diff'); 
		}
		return('more'); # rg1 has more ELEs than rg2; 
	} elsif (@$rg1 < @$rg2) {
		my %h2 = map { $_=>1 } @$rg2; 
		for (@$rg1) {
			defined $h2{$_} or return('diff'); 
		}
		return('less'); 
	} else {
		die "Bad :1:\n"; 
	}
	return(); 
}# compare_grpContent() 

sub len_g2c {
	my ($h) = @_; 
	my @hk = grep { $_ ne 'same' } keys %$h; 
	if      ( scalar(@hk) > 1 ) {
		return('diff'); 
	} elsif ( scalar(@hk) == 1 ) {
		return($hk[0]); 
	} elsif ( @hk == 0 ) {
		return('same'); 
	} else {
		die "Bad :5:\n"; 
	}
}# len_g2c() 


sub type_g2c {
	my ($h) = @_; 
	my @hk = grep { $_ ne 'same' } keys %$h; 
	if      ( scalar(@hk) > 1 ) {
		return('diff'); 
	} elsif ( scalar(@hk) == 1 ) {
		return($hk[0]); 
	} elsif ( @hk == 0 ) {
		return('same'); 
	} else {
		die "Bad :2:\n"; 
	}
	#if      ( defined $h->{'diff'}) {
	#	return('diff'); 
	#} elsif ( scalar(@hk) == 1) {
	#	return( $hk[0] ); 
	#} elsif ( defined $h->{'less'} and defined $h->{'more'} ) {
	#	return('diff'); 
	#} elsif ( scalar(@hk) == 2 ) {
	#	defined $h->{'less'} and defined $h->{'same'} and return('less'); 
	#	defined $h->{'more'} and defined $h->{'same'} and return('more'); 
	#} else {
	#	return('diff'); 
	#}
}# type_g2c()

sub update_grp {
	# If $inGrp, update $grp and return(1); Or else, return(0); 
	my ($grp, $cur, $tag) = @_; 
	$tag //= ''; 
	defined $cur->{'key'} or $cur->{'key'} = [ grep { $_ ne $n_geno } sort keys %{$cur->{'grp_indv'}} ]; 
	my @cur_gType = grep { $_ ne $n_add } @{$cur->{'key'}}; 
	my @grp_gKeys = grep { $_ ne $n_add } @{$grp->[3]}; # Without $n_geno

	my @cur_perms; 
	my $grp_len = scalar(@grp_gKeys); 
	my $cur_len = scalar(@cur_gType); 
	# Make group length no more than current length; 
	if ($grp_len >  $cur_len) {
		for (my $i=0; $i<@grp_gKeys; $i++) {
			$cur_gType[$i] //= $n_add; 
			$cur->{'grp_indv'}{ $cur_gType[$i] } //= [];  
		}
		$cur_len = scalar(@cur_gType); 
	}
	# Compare groups by group length; 
	my %compare_g2c_len; # Compare grp and cur by ELE number in groups; 
	my (%good_grp); 
	for (my $i=0; $i<@grp_gKeys; $i++) {
		for (my $j=0; $j<@cur_gType; $j++) {
			$compare_g2c_len{$i}{$j} = &compare_grpLen( $grp->[2]{$grp_gKeys[$i]}, $cur->{'grp_indv'}{$cur_gType[$j]} ); 
		}
	}
	@cur_perms = &mathSunhh::permutations([0 .. $#cur_gType], $grp_len); 
	for my $cur_perm (@cur_perms) {
		my %len_type; 
		for (my $i=0; $i<$grp_len; $i++) {
			my $j=$cur_perm->[$i]; 
			$len_type{ $compare_g2c_len{$i}{$j} } = 1; 
		}
		my $this_g2c_len = &len_g2c(\%len_type); 

		$this_g2c_len eq 'diff' and next; 
		$grp_len < $cur_len  and $this_g2c_len eq 'more' and next; 
		my $k = join(':', @$cur_perm); $good_grp{$k} = 1; 
	}
	scalar(keys %good_grp) == 0 and return(0); 

	# Compare groups by group elements; 
	my %compare_g2c_ele; # This is required for every grp-cur pair. 
	# For saving time, I'll do this comparison when needed. 
	for my $t1 (keys %good_grp) {
		my %ele_type; 
		my @cur_idx = split(/:/, $t1); 
		for (my $i=0; $i<@cur_idx; $i++) {
			# $i stands for grp_idx; 
			# $j stands for cur_idx; 
			my $j=$cur_idx[$i]; 
			$compare_g2c_ele{$i}{$j} //= &compare_grpContent( $grp->[2]{$grp_gKeys[$i]}, $cur->{'grp_indv'}{$cur_gType[$j]} ); 
			$ele_type{ $compare_g2c_ele{$i}{$j} } = 1; 
		}
		my $this_g2c_ele = &type_g2c(\%ele_type); # Should be one of diff|same|more|less; 
		$this_g2c_ele eq 'diff' and next; 
		if      ($grp_len == $cur_len) {
			# Stored group has equal or more ELEs than current line, so we accept 'same|more|less'; 
			push(@{$grp->[0]}, @{$cur->{'loc'}}); 
			if      ($this_g2c_ele eq 'same' or $this_g2c_ele eq 'more') {
				# Stored group is no less than current line, so we want to add current to stored; 
				# Nothing to do. 
			} elsif ($this_g2c_ele eq 'less') {
				# Stored group is less than current line, so we want to update stored group; 
				$grp->[1] = [@{$cur->{'line'}}]; 
				$grp->[2] = { %{$cur->{'grp_indv'}} }; 
				$grp->[3] = [ grep { $_ ne $n_add } @cur_gType ]; 
				$grp->[4] = [ @{$cur->{'present'}} ]; 
			} else {
				die "Bad :3:\n"; 
			}
			return(1); 
		} elsif ($grp_len < $cur_len) {
			# Stored group has less groups than current line, so only grp <= cur is accepted. 
			$this_g2c_ele eq 'more' and next; 
			push(@{$grp->[0]}, @{$cur->{'loc'}}); 
			$grp->[1] = [@{$cur->{'line'}}]; 
			$grp->[2] = { %{$cur->{'grp_indv'}} }; 
			$grp->[3] = [ grep { $_ ne $n_add } @cur_gType ]; 
			$grp->[4] = [ @{$cur->{'present'}} ]; 
			return(1); 
		} else {
			die "Bad :4:grp_len=$grp_len;cur_len=$cur_len\n"; 
		}
	}# for my $t1 (keys %good_grp) 
	return(0); 
}# update_grp() 


