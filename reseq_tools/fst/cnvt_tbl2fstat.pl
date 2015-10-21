#!/usr/bin/perl 
# Convert snp.tbl format to hierfstat input table, in which N should be writen as NA. 
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use SNP_tbl; 
my $st_obj = SNP_tbl->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"startColN:i", # 2 
	"noHeader!", 
	"ind2grp_list:s", "indCol!", 
	"o_mrk_info:s", # mrk_info
); 

$opts{'startColN'} //= 2; 
$opts{'o_mrk_info'} //= 'mrk_info'; 

my $help_txt = <<HH; 

perl $0 in.snp -ind2grp_list grp_lsit > in.snp.fstat_table

-help 
-ind2grp_list     [list_indv_to_group]Required. Format : IndvID \\t GrpID
-indCol           [Boolean] The IndvID is col_number if given. 
-o_mrk_info       [$opts{'o_mrk_info'}]
-startColN        [$opts{'startColN'}] 

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my %ind2grp = %{ &load_ind2grp( $opts{'ind2grp_list'} ) }; 
my @use_col ; 
my @grpIDs; 
my %allele2num = qw(A 1 C 2 G 3 T 4 N NA); 

sub load_ind2grp {
	my %h; 
	my $fh = &openFH( $_[0], '<' ); 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		$h{$ta[0]} = $ta[1]; 
	}
	close($fh); 
	return \%h; 
}# load_ind2grp () 

unless ($opts{'noHeader'}) {
	my $header = <>; 
	chomp($header); 
	my @ta = split(/\t/, $header); 
	for (my $i=0; $i<@ta; $i++) {
		defined $ind2grp{$ta[$i]} or next; 
		push( @use_col, $i ); 
		push( @grpIDs, $ind2grp{$ta[$i]} ); 
	}
} elsif ($opts{'indCol'}) {
	for ( sort {$a<=>$b} keys %ind2grp ) {
		push( @use_col, $_ ); 
		push( @grpIDs, $ind2grp{$_} ); 
	}
} else {
	&stopErr("[Err] No header or indCol assigned.\n"); 
}
scalar(@use_col) > 0 or &stopErr("[Err] There is no column to be used.\n"); 

my @geno; 
# [Sunhh@Penguin test]$ more ../tt.snp.fmt.1
# grpID    L1      L2      L3      L4      L5      L6      L7      L8      L9
# 1        11      11      22      22      22      22      33      11      22


my $om_fh = &openFH($opts{'o_mrk_info'}, '>'); 
my @mrkIDs; 
my $mrk_copy; 
my $mrk_num = 0; 
while (<>) {
	$. % 10e3 == 1 and &tsmsg("[Msg] $. line.\n"); 
	chomp; 
	my @ta = split(/\t/, $_); 
	my $tb = &geno2num( [ @ta[@use_col] ] ); 
	$mrk_num ++; 
	push(@mrkIDs, "L$mrk_num"); 
	print {$om_fh} "$mrkIDs[-1]\t$ta[0]\t$ta[1]\n"; 
	$mrk_copy //= "$mrkIDs[-1]r\t$ta[0]\t$ta[1]\n"; 
	for (my $i=0; $i<@grpIDs; $i++) {
		push(@{$geno[$i]}, $tb->[$i]); 
	}
}
if ($mrk_num == 1) {
	&tsmsg("[Wrn] Only 1 site in the data, with is not enough.\n"); 
	print {$om_fh} $mrk_copy; 
	push(@mrkIDs, "$mrkIDs[0]r"); 
}
if ($mrk_num == 0) {
	&stopErr("[Err] There isn't any available site.\n"); 
}
close($om_fh); 

print STDOUT join("\t", 'GrpID', @mrkIDs)."\n"; 

for (my $i=0; $i<@grpIDs; $i++) {
	if ( $mrk_num == 1 ) {
		@{$geno[$i]} == 1 or &stopErr("[Err] geno=[@{$geno[$i]}], but mrk_num == 1\n"); 
		@{$geno[$i]} = (@{$geno[$i]}, @{$geno[$i]}); # repeat once. 
	}
	print join("\t", $grpIDs[$i], @{$geno[$i]})."\n"; 
}
&tsmsg("[Rec] done.\n"); 

sub geno2num {
	# $_[0] : [$geno1, $geno2, $geno3, ... ]
	# $_[1] : Ploidy number for MAF counting. 2 by default. Could be only 1 or 2. 
	# Only count ploidy == 1 or 2 ; 
	$_[1] //= 2; 
	$_[1] == 2 or $_[1] == 1 or &stopErr("[Err] Bad ploidy number [$_[1]] which could be only 1 or 2.\n"); 
	my %cnt; 
	my @new_bb; 
	for (@{$_[0]}) {
		($_ eq 'N' or $_ eq 'n') and do { $cnt{'N'}++; push(@new_bb, ['N', 'N']); next; }; 
		if ( $_ =~ m/^[ATGC*]$/i or $_ =~ m/\+/ ) {
			$cnt{'allele2cnt'}{$_} += $_[1]; 
			push(@new_bb, [$_, $_]); 
		} elsif ( $_ =~ m/^([ATGC])\*$/i or $_ =~ m/^\*([ATGC])$/i ) {
			$_[1] == 1 and do { &tsmsg("[Wrn] Bad genotype [$_] for ploidy=1.\n"); $cnt{'N'}++; push(@new_bb, ['N', 'N']); next; }; 
			$cnt{'allele2cnt'}{ uc($1) } ++; 
			$cnt{'allele2cnt'}{'*'} ++; 
			push(@new_bb, [uc($1), '*']); 
		} elsif ( ( my @bb = &SNP_tbl::dna_d2b( &SNP_tbl::dna_b2d($_) ) ) > 0 ) {
			@bb   == 1 and do { $cnt{'allele2cnt'}{$bb[0]} += $_[1]; push(@new_bb, [$bb[0], $bb[0]]); next; }; 
			$_[1] == 1 and do { &tsmsg("[Wrn] Bad genotype [$_] for ploidy=1.\n"); $cnt{'N'}++; push(@new_bb, ['N', 'N']); next; }; 
			@bb   >  2 and do { &tsmsg("[Wrn] Bad genotype [$_] for ploidy=2.\n"); $cnt{'N'}++; push(@new_bb, ['N', 'N']); next; }; 
			for my $tb ( @bb ) { $cnt{'allele2cnt'}{$tb}++; } 
			push(@new_bb, [@bb]); 
		} else {
			&tsmsg("[Wrn] Weired genotype [$_] is treated as homozygous.\n"); 
			$cnt{'allele2cnt'}{$_} += $_[1]; 
			push(@new_bb, [$_, $_]); 
		}
	}# End for (@{$_[0]})

	$cnt{'N'} //= 0; 
	my @aa = sort { $cnt{'allele2cnt'}{$b} <=> $cnt{'allele2cnt'}{$a} || $a cmp $b } keys %{$cnt{'allele2cnt'}}; 
	my $nextN = 5; 
	my %a2n = %allele2num; 
	for (@aa) {
		defined $a2n{$_} and next; 
		$a2n{$_} = ( $nextN < 10 ) ? $nextN : 'NA' ; 
		$nextN++; 
	}
	my @back_num; 
	for (@new_bb) {
		( $a2n{ $_->[0] } eq 'NA' or $a2n{ $_->[1] } eq 'NA' ) and do { push(@back_num, 'NA'); next; }; 
		push( @back_num, join('', $a2n{ $_->[0] }, $a2n{ $_->[1] }) ); 
	}

	return (\@back_num); 
}# sub geno2num() 



