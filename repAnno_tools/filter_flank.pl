#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 

# dgt.tab.wPP.filtN 
# eleID   eleS    eleE    Str     seqID   LTR1_S  LTR1_E  LTR2_S  LTR2_E  Inner_S Inner_E PBS_S   PBS_E   PPT_S   PPT_E   scfID
# RR2     785204  788754  +       seq10   785209  785369  788589  788749  785370  788588  785373  785384  -1      -1      S400016_pilon
# RR3     839740  842441  +       seq10   839745  840173  842011  842436  840174  842010  840174  840184  -1      -1      S400016_pilon

# raw_scf.fasta

!@ARGV and die "perl $0 srcFa dgt.tab.wPP.filtN\n"; 

my $srcFa = shift; 
my $exeMuscle = "/data/Sunhh/src/Align/muscle/muscle3.8.31_i86linux64"; 

my $flank_len = 50; 
my $min_ident_pos = int($flank_len/2); 
my $min_ident_rat = 0.6; 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ( $ta[0] eq 'eleID' ) {
		print STDOUT join("\t", $_, qw/LTR1_cov LTR1_ident LTR2_cov LTR2_ident/)."\n"; 
		next; 
	}
	my ($scfID, $ltr1S, $ltr1E, $ltr2S, $ltr2E) = @ta[ 15, 5,6, 7,8 ]; 
	# For 5p upstream
	my %cnt1 = %{&aln_flank( $scfID, $ltr1S-$flank_len, $ltr1S-1, $scfID, $ltr2S-$flank_len, $ltr2S-1 )}; 
	# For 3p downstream
	my %cnt2 = %{&aln_flank( $scfID, $ltr1E+1, $ltr1E+$flank_len, $scfID, $ltr2E+1, $ltr2E+$flank_len )}; 
	print STDOUT join("\t", $_, $cnt1{id_pos}/$cnt1{l1}, $cnt1{identity}, $cnt2{id_pos}/$cnt2{l1}, $cnt2{identity})."\n"; 
}

# print join("\t", $_, $cnt{id_pos}, $cnt{id_pos}-$cnt{idN_pos}, $cnt{aln_pos}, $cnt{identity})."\n"; 
sub aln_flank {
	my ( $id1, $s1, $e1, $id2, $s2, $e2 ) = @_; 
	$s1 < 0 and $s1 = 1; 
	$s2 < 0 and $s2 = 1; 
	&exeCmd("rm -f tmp_in.fa tmp_in.afa"); 
	open UP,'>',"tmp_up.lis" or die; 
	print UP join("\t", $id1, $s1, $e1, "+", "seqUP1LTR")."\n"; 
	print UP join("\t", $id2, $s2, $e2, "+", "seqUP2LTR")."\n"; 
	close UP; 
	&exeCmd("deal_fasta.pl -drawByList -drawLcol 0,1,2,3,4 -drawList tmp_up.lis $srcFa > tmp_in.fa"); 
	&exeCmd("$exeMuscle -in tmp_in.fa -out tmp_aln.afa"); 
	open F,'<',"tmp_aln.afa" or die; 
	my (@id, @seq); 
	while ( my $l=<F> ) {
		chomp; 
		if ( $l =~ m/>(\S+)/ ) {
			push(@id, $1); 
			push(@seq, ""); 
		} else {
			$seq[-1] .= $l; 
		}
	}
	close F; 
	for my $t (@seq) {
		$t =~ s/\s//g; 
		$t = uc($t); 
	}
	scalar(@seq) == 2 and length($seq[0]) == length($seq[1]) or die "fadsf\n$_\n"; 
	my @b1 = split(//, $seq[0]); 
	my @b2 = split(//, $seq[1]); 
	my %cnt; 
	for (my $i=0; $i<@b1; $i++) {
		$b1[$i] eq '-' and next; 
		$cnt{l1}++; 
		$b2[$i] eq '-' and next; 
		$b1[$i] ne 'N' and $b2[$i] ne 'N' and $cnt{aln_pos} ++; 
		$b1[$i] eq $b2[$i] and $cnt{id_pos} ++; 
		$b1[$i] eq 'N' and $b1[$i] eq $b2[$i] and $cnt{idN_pos} ++; 
	}
	for my $t (qw/id_pos idN_pos aln_pos/) {
		defined $cnt{$t} or $cnt{$t} = 0; 
	}
	$cnt{identity} = ($cnt{aln_pos} == 0) ? -1 : ( $cnt{id_pos}-$cnt{idN_pos} )/$cnt{aln_pos}; 
	return \%cnt; 
}
