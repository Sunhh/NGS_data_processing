#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use fastaSunhh; 
my $fas_obj = fastaSunhh->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"agpFile:s", 
	"ctgFas:s", 
	"noAddRest!", 
	"Nchar:s", "Uchar:s", 
); 
$opts{'Nchar'} //= 'n'; 
$opts{'Uchar'} //= 'n'; 

my $help_txt = <<HH; 

perl $0 -agpFile SP_Nov_11_13_30_8_90_3_superscaffold_merged_edit.agp -ctgFas Spo_V1p5.scaf.fa > Spo_V1p5_Jn.scaf.fa

-help 

-noAddRest          [Bool] Only output sequences within .agp file if given. 
-Nchar              ['$opts{'Nchar'}'] Character to fill gap for ele_N. 
-Uchar              ['$opts{'Uchar'}'] Character to fill gap for ele_U. 

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
for my $tk (qw/agpFile ctgFas/) {
	defined $opts{$tk} or &LogInforSunhh::usage($help_txt); 
}
my %orient = qw(
+     P
-     M
?     P
0     P
na    P
P     P
M     M
1     P
-1    M
); 

&tsmsg("[Msg] Loading AGP file [$opts{'agpFile'}].\n"); 
my @agp_arr = @{ &agp_sort( &load_agp_file( $opts{'agpFile'} ) ) }; 
&tsmsg("[Msg] Loading ctg_fasta [$opts{'ctgFas'}].\n"); 
my %ctg_fas = %{ $fas_obj->save_seq_to_hash('faFile'=>$opts{'ctgFas'}, 'has_head'=>1) }; 
for (keys %ctg_fas) { $ctg_fas{$_}{'seq'} =~ s/\s//g; } 

&tsmsg("[Msg] Joining sequences.\n"); 
my %new_seq; 
my $seq_cnt = 0; 
my %used_ctg; 
for my $la (@agp_arr) {
	my $newID = $la->[0]; 
	defined $new_seq{$newID} or $new_seq{$newID}{'order'} = ++$seq_cnt; 
	$new_seq{$newID}{'prevP'} //= 0; 
	$new_seq{$newID}{'seq'} //= ''; 
	my $add_seq = ''; 
	$la->[1] == $new_seq{$newID}{'prevP'} + 1 or &stopErr("[Err] The position coordinates are not continuous: scfID=$newID prevP=$new_seq{$newID}{'prevP'} currP=$la->[1]\n"); 
	if ( $la->[4] =~ m/^N$/i ) {
		$add_seq = $opts{'Nchar'} x $la->[5]; 
		$new_seq{$newID}{'seq'} .= $add_seq; 
	} elsif ( $la->[4] =~ m/^U$/i ) { 
		$add_seq = $opts{'Uchar'} x $la->[5]; 
		$new_seq{$newID}{'seq'} .= $add_seq; 
	} elsif ( $la->[4] eq 'W') {
		defined $ctg_fas{$la->[5]} or &stopErr("No ctg_fas for ID [$la->[5]]\n"); 
		$add_seq = substr( $ctg_fas{$la->[5]}{'seq'}, $la->[6]-1, $la->[7]-$la->[6]+1 ); 
		$orient{ $la->[8] } eq 'M' and &fastaSunhh::rcSeq(\$add_seq, 'rc'); 
		$new_seq{$newID}{'seq'} .= $add_seq; 
		$used_ctg{$la->[5]} = 1; 
	} else {
		&stopErr("[Err] No rule for component_type[5] [$la->[4]] yet.\n"); 
	}
	$new_seq{$newID}{'prevP'} = $la->[2]; 
}

&tsmsg("[Msg] Output joined sequences.\n"); 
for my $tk (sort {$new_seq{$a}{'order'} <=> $new_seq{$b}{'order'}} keys %new_seq) {
	$new_seq{$tk}{'seq'} =~ s/(.{60})/$1\n/g; chomp( $new_seq{$tk}{'seq'} ); 
	print STDOUT ">$tk\n$new_seq{$tk}{'seq'}\n"; 
}

unless ($opts{'noAddRest'}) {
	&tsmsg("[Msg] Output rest sequences.\n"); 
	for my $tk (sort { $ctg_fas{$a}{'Order'} <=> $ctg_fas{$b}{'Order'} } keys %ctg_fas) {
		defined $used_ctg{$tk} and next; 
		$ctg_fas{$tk}{'seq'} =~ s/(.{60})/$1\n/g; chomp( $ctg_fas{$tk}{'seq'} ); 
		print STDOUT ">$tk\n$ctg_fas{$tk}{'seq'}\n"; 
	}
}

&tsmsg("[Msg] Done $0.\n"); 

# [Sunhh@whale super_scaffold]$ head -3 SP_Nov_11_13_30_8_90_3_superscaffold_merged.agp
# ##agp-version   2.0
# Super_scaffold_248      1       341379  1       W       SpoScf_00639    1       341379  -
# Super_scaffold_248      341380  600696  2       N       259317  scaffold        yes     map
sub load_agp_file {
	# https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/ 
	my $fh = &openFH($_[0], '<'); 
	my @back; 
	while (<$fh>) {
		m/^\s*(#|$)/ and next; 
		chomp; 
		my @ta = split(/\t/, $_); 
		push(@back, [@ta]); 
	}
	close($fh); 
	return(\@back); 
}

sub agp_sort {
	# $_[0] : \@agp_arr 
	my @back; 
	my %ord; 
	my $cnt=0; 
	for (@{$_[0]}) {
		$ord{$_->[0]} //= ++$cnt; 
	}
	@back = sort { $ord{$a->[0]} <=> $ord{$b->[0]} || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @{$_[0]}; 
	return(\@back); 
}# agp_sort() 

