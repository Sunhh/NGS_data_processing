#!/usr/bin/perl
use strict; 
use warnings; 

use LogInforSunhh; 
use fileSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"debug!", 
	"includeEle:s@", "chkEle:s@", 
	"sepEle:s@", 
	"outDupPair:s", 
	"dupSrcGff:s", 
); 

sub usage {
	print <<HH; 
################################################################################
# perl $0 dup_aug1_aug2.gff3 > nonDup_aug1_aug2.gff3
#
# -help 
# -debug
#
# -includeEle  [match match_part]
# -sepEle      [match]
# -chkEle      [Same to includeEle]
#
# -outDupPair  [Null] Output file recording duplicated pair of gene models. 
#
# -dupSrcGff  [Null] Check duplication from this file if given. 
#
################################################################################
HH
	exit 1; 
}

-t and !@ARGV and &usage(); 
$opts{'help'} and &usage(); 
sub mk_ele {
	my ($hr, $ar1, $ar2) = @_; 
	my $aa = 0; 
	if ( defined $ar1 and scalar(@$ar1) > 0 ) {
		%$hr = map { lc($_) => $aa++ } @$ar1; 
	} elsif ( defined $ar2 and scalar(@$ar2) > 0 ) {
		%$hr = map { lc($_) => $aa++ } @$ar2; 
	} else {
		; 
	}
}

my (%ele2incl, %ele2stop, %ele2chk); 
&mk_ele( \%ele2incl, $opts{'includeEle'}, [ qw/match match_part/ ] ); 
&mk_ele( \%ele2stop, $opts{'sepEle'},     [ qw/match/ ]); 
&mk_ele( \%ele2chk,  $opts{'chkEle'},     [ sort { $ele2incl{$a} <=> $ele2incl{$b} } keys %ele2incl ]); 

my $oDupPairFh; 
if ( defined $opts{'outDupPair'} ) {
	open $oDupPairFh, '>', "$opts{'outDupPair'}" or &stopErr("[Err] Unavailable to write file [$opts{'outDupPair'}]\n"); 
}

my %have_model; 

# ##gff-version 3
# S401083_pilon   AUGUSTUS        match   113     6130    0.39    -       .       ID=0:g1.t1;
# S401083_pilon   AUGUSTUS        match_part      113     472     1       -       0       ID=0:g1.t1.c1;Parent=0:g1.t1;
# S401083_pilon   AUGUSTUS        match_part      595     807     0.99    -       0       ID=0:g1.t1.c2;Parent=0:g1.t1;
# S401083_pilon   AUGUSTUS        match_part      1269    1412    0.91    -       0       ID=0:g1.t1.c3;Parent=0:g1.t1;

if (defined $opts{'dupSrcGff'}) {
	my $dFh = &openFH($opts{'dupSrcGff'}, '<'); 
	my @geneLines; 
	while (<$dFh>) {
		m/^#/ and next; 
		m/^\s*$/ and next; 
		m/^>/ and last; 
		chomp; 
		my @ta = split(/\t/, $_); 
		my $eleType = lc($ta[2]); 
		defined $ele2incl{$eleType} or next; 
		if ( defined $ele2stop{$eleType} ) { 
			if ( scalar(@geneLines) > 0 ) {
				&add_mod(\@geneLines, \%have_model); 
			}
			@geneLines = (); 
			push(@geneLines, [[@ta]]); 
		} elsif ( defined $ele2incl{$eleType} ) {
			push( @geneLines, [[@ta]] ); 
		} else {
			die; 
		}
	}
	close $dFh; 
	scalar(@geneLines) > 0 and &add_mod(\@geneLines, \%have_model); 
}

my @geneLines; 
while (<>) {
	m/^#/ and next; 
	m/^\s*$/ and next; 
	m/^>/ and last; 
	chomp; 
	my @ta = split(/\t/, $_); 
	my $eleType = lc( $ta[2] ); 
	defined $ele2incl{$eleType} or next; 
	if ( defined $ele2stop{$eleType} ) {
		if (scalar(@geneLines) > 0) {
			my $noAdd = (defined $opts{'dupSrcGff'}) ? 1 : 0 ; 
			&is_have(\@geneLines, \%have_model, $noAdd) or &out_gL(\@geneLines); 
		}
		@geneLines = (); 
		push(@geneLines, [[@ta]]); 
	} elsif ( defined $ele2incl{$eleType} ) {
		push(@geneLines, [[@ta]]); 
	} else {
		die; 
	}
}

if (scalar(@geneLines) > 0) {
	my $noAdd = (defined $opts{'dupSrcGff'}) ? 1 : 0 ;
	&is_have(\@geneLines, \%have_model, $noAdd) or &out_gL(\@geneLines); 
	@geneLines = (); 
}

sub add_mod {
	my ($ar, $hr) = @_; 
	my $scfID = $ar->[0][0][0];
	my ($gS, $gE, $gStr) = @{$ar->[0][0]}[3,4,6]; 
	my $gKey ; 
	my @chk_ar ;
	my ($chkS, $chkE);
	for my $t1 (@$ar) {
		my $eleType = lc($t1->[0][2]);
		defined $ele2chk{$eleType} or next;
		push(@chk_ar, $t1);
		defined $chkS or $chkS = $t1->[0][3];
		defined $chkE or $chkE = $t1->[0][4];
		$chkS > $t1->[0][3] and $chkS = $t1->[0][3];
		$chkE < $t1->[0][4] and $chkE = $t1->[0][4];
	}
	$gKey = join("\t", $chkS, $chkE, $gStr);
	push(@{$hr->{$scfID}{$gKey}}, [@$ar]); 
}

sub is_have {
	my ($ar, $hr) = @_; 
	my $noAdd = shift; 
	$noAdd = $noAdd // 0; 
	my $scfID = $ar->[0][0][0]; 
	my $is_have = 1; 
	my ($gS, $gE, $gStr) = @{$ar->[0][0]}[3,4,6]; 
	my $gKey = join("\t", $gS, $gE, $gStr); 
	my $had_ta8 = ''; 

	my @chk_ar ; 
	my ($chkS, $chkE); 
	for my $t1 (@$ar) {
		my $eleType = lc($t1->[0][2]); 
		defined $ele2chk{$eleType} or next; 
		push(@chk_ar, $t1); 
		defined $chkS or $chkS = $t1->[0][3]; 
		defined $chkE or $chkE = $t1->[0][4]; 
		$chkS > $t1->[0][3] and $chkS = $t1->[0][3]; 
		$chkE < $t1->[0][4] and $chkE = $t1->[0][4]; 
	}
	# Update gKey
	$gKey = join("\t", $chkS, $chkE, $gStr); 
	
	if ( defined $hr->{$scfID} ) {
		if ( defined $hr->{$scfID}{$gKey} ) {
			$is_have = 0; 
			CHK_MODEL: 
			for my $tr (@{$hr->{$scfID}{$gKey}}) {
				my $tmp_have = 1; 
				my @chk_tr; 
				for my $t1 (@$tr) {
					my $eleType = lc($t1->[0][2]); 
					defined $ele2chk{$eleType} or next; 
					push(@chk_tr, $t1); 
				}
				scalar( @chk_ar ) == scalar( @chk_tr ) or do { $tmp_have = 0; next CHK_MODEL; }; 
				CHK_EXON: 
				for (my $i=0; $i<@chk_tr; $i++) {
					$chk_ar[$i][0][3] eq $chk_tr[$i][0][3] or do { $tmp_have = 0; last CHK_EXON; } ; 
					$chk_ar[$i][0][4] eq $chk_tr[$i][0][4] or do { $tmp_have = 0; last CHK_EXON; } ; 
					$chk_ar[$i][0][6] eq $chk_tr[$i][0][6] or do { $tmp_have = 0; last CHK_EXON; } ; 
				}
				if ( $tmp_have == 1 ) {
					$is_have = 1; 
					$had_ta8 = $tr->[0][0][8]; 
					defined $opts{'outDupPair'} and &out_gL( $tr, $oDupPairFh ); 
					defined $opts{'outDupPair'} and &out_gL( $ar, $oDupPairFh ); 
					last CHK_MODEL; 
				}
			}
		} else {
			$is_have = 0; 
		}
	} else {
		$is_have = 0; 
	}
	if ($is_have == 0) {
		$noAdd or push(@{$hr->{$scfID}{$gKey}}, [@$ar]); 
	} else {
		$opts{'debug'} and &tsmsg("[Msg] Skip duplicted gKey [$gKey] for [$ar->[0][0][8]] by [$had_ta8]\n"); 
	}
	
	return $is_have; 
}

sub out_gL {
	my $ar = shift; 
	my $fh = shift; 
	defined $fh or $fh = \*STDOUT; 
	for my $r1 (@$ar) {
		print {$fh} join("\t", @{$r1->[0]})."\n";
	}
	return 0;
}

