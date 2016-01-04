#!/usr/bin/perl -w 
# 2013-03-05 Add -rawID option. 
use strict; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
		"help!", 
		"frame:s", # Edit here. 
		"seqLen:i", 
		"trimStop!", 
		"noStopCodon!", 
		"rawID!", 
);
my $info = <<INFO;
perl $0 cds
-frame         1,2,3,-1,-2,-3
-seqLen        60
-trimStop      for AA   output sequences.
-noStopCodon   for Nucl output sequences.
-rawID         Use same ID with the CDS sequences.
INFO
($opts{help} or !@ARGV) and die $info; 
(defined $opts{seqLen} and int($opts{seqLen}) > 0) or $opts{seqLen} = 60; 
our %good_sframe; 
%good_sframe = qw(
		1 1 
		2 1
		3 1 
	 -1 -1 
	 -2 -1
	 -3 -1
); 
my @ffs = (1 .. 3); 
if (defined $opts{frame}) {
	@ffs = split(/,/, $opts{frame});
	$opts{frame} == 6 and @ffs = qw(1 2 3 -1 -2 -3); 
	my %used; 
	my @ffs_save; 
	for (@ffs) {
		s/\s//g; 
		defined $good_sframe{$_} or do { warn "Skip unknown sframe [$_].\n"; next; }; 
		defined $used{$_} or push(@ffs_save, $_); 
		$used{$_} = 1; 
	}
	@ffs = @ffs_save; 
}

my @aas = split(//, 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG-'); 
my @sss = split(//, '---M---------------M---------------M-----------------------------'); 
my @bb1 = split(//, 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG-'); 
my @bb2 = split(//, 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG-'); 
my @bb3 = split(//, 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG-'); 

my %cod; 
my %ini; 
my %head; 
{
 for (my $i=0; $i<@aas; $i++) {
  my $bbb = "$bb1[$i]$bb2[$i]$bb3[$i]"; 
  $bbb = uc($bbb); 
  $aas[$i] = uc($aas[$i]); 
  $cod{$bbb} = $aas[$i]; 
  $sss[$i] eq 'M' and $ini{$bbb} = 1; 
 }
}

my (%seq, $key, @kks); 
while (<>) {
 if (s/^>(\S+)//) {
  $key = $1; 
  chomp; 
  $head{$key} = $_; 
  push(@kks, $key); 
 }else {
  $seq{$key} .= $_; 
 }
}

#for my $sframe ( 0..2 ) { 
for my $sframe ( @ffs ) { 
for (@kks) {
 $seq{$_} =~ s/\s//g; 
 $seq{$_} = uc($seq{$_}); 
 my $l = length($seq{$_}); 
 my $tail_tag = ($opts{rawID}) ? '' : ":$sframe" ; 
 print STDOUT ">$_${tail_tag} [sframe=$sframe]$head{$_}\n"; 
 print STDERR ">$_${tail_tag} [sframe=$sframe]$head{$_}\n"; 
 my $aa = ''; 
 my $nn = ''; 
 if ($sframe > 0) {
  for (my $i=$sframe-1; $i<$l; $i+=3) {
   my $e = $i+3; $e > $l and last; 
   my $bbb = substr($seq{$_}, $i, 3); 
   $aa .= (defined $cod{$bbb}) ? $cod{$bbb} : 'X' ; 
   $nn .= $bbb; 
  }
 }elsif ($sframe < 0) {
	my $rc_seq = reverse($seq{$_}); 
	$rc_seq =~ tr/acgturykmbvdhACGTURYKMBVDH/tgcaayrmkvbhdTGCAAYRMKVBHD/; 
	for (my $i=-$sframe-1; $i<$l; $i+=3) {
		my $e = $i+3; $e > $l and last; 
		my $bbb = substr($rc_seq, $i, 3); 
		$aa .= (defined $cod{$bbb}) ? $cod{$bbb} : 'X'; 
		$nn .= $bbb; 
	}
 }else {
	die "sframe should not be 0.\n"; 
 }
 $opts{trimStop} and $aa =~ s/\*$//; 
 if ( $opts{noStopCodon} ) {
	my $tb = substr($nn, -3, 3); 
	$cod{$tb} eq '*' and substr($nn, -3, 3) = ''; 
 }
 if ($opts{seqLen} < 32766) {
  $aa =~ s/(.{$opts{seqLen}})/$1\n/og; 
  chomp($aa); 
  print STDOUT "$aa\n"; 
  $nn =~ s/(.{$opts{seqLen}})/$1\n/go;
  chomp($nn); 
  print STDERR "$nn\n"; 
 }else{
	my $aa_l = length($aa); 
	my $nn_l = length($nn); 
	for (my $j=0; $j<$aa_l; $j+=$opts{seqLen}) {
		my $e = $j+$opts{seqLen}; 
		$e > $aa_l and $e=$aa_l; 
		my $sub_oseq = substr($aa, $j, $e-$j); 
		print STDOUT "$sub_oseq\n"; 
	}
	for (my $j=0; $j<$nn_l; $j+=$opts{seqLen}) {
		my $e = $j+$opts{seqLen}; 
		$e > $nn_l and $e=$nn_l; 
		my $sub_oseq = substr($nn, $j, $e-$j); 
		print STDERR "$sub_oseq\n"; 
	}
 }
}
}




