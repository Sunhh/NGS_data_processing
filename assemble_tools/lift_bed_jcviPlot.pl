#!/usr/bin/perl
# [Sunhh@panda aa]$ ColLink.pl Ma.bed -f1 Ma.lifted.bed -keyC1 3,4 -keyC2 3,4 -add -Col1 0-2 2>2 | head -4
# Cma_Scf00001    25078   25079   pop2-4:1.653495 Cma_Scf00001:25079      chr4    25078   25079
# Cma_Scf00001    25079   25080   pop2-4:1.653495 Cma_Scf00001:25080      chr4    25079   25080
# Cma_Scf00001    41840   41841   pop2-4:0.000000 Cma_Scf00001:41841      chr4    41840   41841
# [Sunhh@panda aa]$ head -4 Ma.lifted.bed
# chr1    20324   20325   pop2-1:13.326251        Cma_Scf00018:20325
# chr1    24767   24768   pop2-1:13.326251        Cma_Scf00018:24768
use strict; 
use warnings; 
use LogInforSunhh; 
use mathSunhh; 
my $ms_obj = mathSunhh->new();
use fileSunhh; 

my $help_txt = <<HH; 

perl $0 ABC.agp ABC.bed | deal_table.pl -column 0-4 > ABC.lifted.bed

# Written for converting ALLMAPs ABC.bed file to ABC.lifted.bed file according to ABC.agp file.
# After preparing this, using files ABC.bed, ABC.lifted.bed and ABC.agp, we can plot chromosomems by command :
#   python -m jcvi.assembly.allmaps plot ABC.bed chromID
# The default output is chromID.pdf .

HH
!@ARGV and &LogInforSunhh::usage($help_txt); 

my $f_agp = shift; 
my %agp_c2s = %{ &fileSunhh::load_agpFile( $f_agp ) }; 

my @out; 
while (<>) {
	m/^\s*(#|$)/ and next; 
	chomp; 
	my @ta = &splitL("\t", $_); 
	my $ctg_ID  = $ta[0]; 
	my $ctg_pos = $ta[2]; 
	$ta[2]-$ta[1] == 1 or die "$_\n"; 
	my $ctg_str = '+'; 
	my @new_scfInf = $ms_obj->transfer_position( 'from_ref2qry' => {}, 'to_qry2ref' => \%agp_c2s, 'fromLoc' => [$ctg_ID, $ctg_pos, $ctg_str] ); 
	my $scf_ID  = $new_scfInf[0]; 
	my $scf_pos = $new_scfInf[1]; 
	push(@out, [@ta, $scf_ID, $scf_pos]); 
}

# for my $tr1 (sort { &srt_chrN($a->[-2], $b->[-2]) || $a->[-1] <=> $b->[-1] } @out) {
# for my $tr1 (sort { $a->[-2] <=> $b->[-2] || $a->[-2] cmp $b->[-2] || $a->[-1] <=> $b->[-1] } @out) {
for my $tr1 (sort { $a->[-2] cmp $b->[-2] || $a->[-1] <=> $b->[-1] } @out) {
	my @nn = (3 .. ($#$tr1-2), 0, 1, 2); 
	print join("\t", $tr1->[$#$tr1-1], $tr1->[$#$tr1]-1, $tr1->[$#$tr1], @{$tr1}[@nn])."\n"; 
}

sub srt_chrN {
	my ($i1, $i2) = @_; 
	my ($n1, $n2); 
	$i1 =~ m/^chr[^\s\d]*(\d+)$/i and $n1 = $1; 
	$i2 =~ m/^chr[^\s\d]*(\d+)$/i and $n2 = $1; 
	$n1 //= $i1; 
	$n2 //= $i2; 
	return( $n1 <=> $n2 || $n1 cmp $n2 ); 
}

