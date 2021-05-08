#!/usr/bin/perl
use strict;
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"opref:s", 
	'inCds:s', 
);
my $help_txt = <<HH; 
################################################################################
# perl $0 -opref outPrefix   -inCds all_cds.fasta input.trnascanse.out
################################################################################
HH
$opts{'help'} and &LogInforSunhh::usage($help_txt); 
-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

$opts{'opref'} //= 'out'; 
my %ofn; 
&init_ofiles(); 


my (%cds, %codonInCDS); 
defined $opts{'inCds'} and &load_cds();

# Sequence                tRNA            Bounds          tRNA    Anti    Intron Bounds   Inf     HMM     2'Str   Hit     Isotype Isotype
# Name            tRNA #  Begin           End             Type    Codon   Begin   End     Score   Score   Score   Origin  CM      Score   Note
# --------        ------  -----           ------          ----    -----   -----   ----    ------  -----   -----   ------  ------- ------- ------
# WM97pbV1_Chr01  1       494557          494628          iMet    CAT     0       0       60.9    31.70   29.20   Inf     iMet    115.6
# WM97pbV1_Chr01  2       816817          816900          Leu     CAA     0       0       79.1    54.80   24.30   Inf     Leu     124.7
# WM97pbV1_Chr01  3       2012293         2012364         Gly     TCC     0       0       73.9    46.20   27.70   Inf     Gly     85.6
### In GtRNAdb, 'tRNA Type' information is used as 'tRNA isotype', so I don't care about the 'Isotype CM' information; 

my %ccc = &iso_codon_table(); 

my %h; 
while (<>) {
	chomp; 
	m!^\S+\s+(\d+)\s+! or next; 
	my @ta = &splitL("\t", $_); 
	my ($isotype, $anticodon) = @ta[4,5]; # isotype is AA type; 
	$h{'antiC'}{$anticodon} ++; 
	$h{'isotype'}{$isotype} ++; 
	$h{'iso2antiC'}{$isotype}{$anticodon} ++; 
}
for my $boxID (qw/FourBox SixBox TwoBox OtherBox/) {
	for my $aa (@{$ccc{$boxID}}) {
		my $ttl_cnt = 0; 
		my %codonV; 
		my $ttl_codonCnt = 0; 
		my $ttl_codonRat = 0; 
		for my $antiC (@{$ccc{'aa2antiC'}{$aa}}) {
			$h{'antiC'}{$antiC} //= 0; 
			$ttl_cnt += $h{'antiC'}{$antiC}; 
			my $c = $ccc{'antiC2C'}{$antiC}; 
			defined $c or die "failed [$antiC]\n"; 
			if ( defined $opts{'inCds'} ) {
				$codonV{$antiC}{'cnt'}   = (defined $codonInCDS{'codonCnt'}{$c})    ? $codonInCDS{'codonCnt'}{$c}    : 0 ; 
				$codonV{$antiC}{'ratio'} = (defined $codonInCDS{'codonRatio'}{$c}) ? $codonInCDS{'codonRatio'}{$c} : 0 ; 
				$ttl_codonCnt += $codonV{$antiC}{'cnt'}; 
				$ttl_codonRat += $codonV{$antiC}{'ratio'}; 
				$codonV{$antiC}{'ratio'} = sprintf("%.2f",  $codonV{$antiC}{'ratio'} * 100); 
			}
		}
		&fileSunhh::write2file($ofn{'boxTbl'}, join("\t", $boxID, $aa, $ttl_cnt, map { "$_ ($h{'antiC'}{$_})" } @{$ccc{'aa2antiC'}{$aa}})."\n", '>>'); 
		if ( defined $opts{'inCds'} ) {
			$ttl_codonRat = sprintf("%.2f", $ttl_codonRat * 100); 
			&fileSunhh::write2file(
				$ofn{'codonTbl'}, 
				join("\t", $boxID, $aa, 
					"$ttl_codonCnt ($ttl_codonRat\%)", 
					( map { "$ccc{'antiC2C'}{$_} ($codonV{$_}{'cnt'} $codonV{$_}{'ratio'}\%)" } @{$ccc{'aa2antiC'}{$aa}} ), 
				)."\n", 
				'>>'
			); 
			for my $ac ( @{$ccc{'aa2antiC'}{$aa}} ) {
				my $c = $ccc{'antiC2C'}{$ac}; 
				&fileSunhh::write2file(
					$ofn{'codon2gene'}, 
					join("\t", 
						$aa, 
						$ac, 
						$c, 
						scalar(keys %{$codonInCDS{'codon2gene'}{$c}}), 
						join(";", sort keys %{$codonInCDS{'codon2gene'}{$c}})
					)."\n", 
					'>>'
				); 
			}
		}
	}
}


# Return : (%back); 
#  %{$back{'aa2antiC'}} = { 'Ala'=>[qw/AGC GGC CGC TGC/], ... }; 
#  @{$back{'FourBox'}}  = [ qw/Ala Gly Pro Thr Val/ ]; 
#  @{$back{'SixBox'}}   = [ qw/Ser Arg Leu/ ]; 
#  @{$back{'TwoBox'}}   = [ qw/Phe Asn Lys Asp Glu His Gln/ ]; 
#  @{$back{'OtherBox'}} = [ qw/Ile Met Tyr Supres Cys Trp SelCys/ ]; 
#  %{$back{'antiC2aa'}} = { 'AGC'=>'Ala', 'GGC'=>'Ala', ... }; 
#  %{$back{'antiC2C'}}  = { 'AGC'=>'GCT', ... }; 
#  %{$back{'C2antiC'}}  = { 'GCT'=>'AGC', ... }; 
sub iso_codon_table {
	my %back; 
	# Edit for Euk : Standard code : transl_table=1 
	@{$back{'aa2antiC'}{'Ala'}} = qw/AGC GGC CGC TGC/; 
	@{$back{'aa2antiC'}{'Gly'}} = qw/ACC GCC CCC TCC/; 
	@{$back{'aa2antiC'}{'Pro'}} = qw/AGG GGG CGG TGG/; 
	@{$back{'aa2antiC'}{'Thr'}} = qw/AGT GGT CGT TGT/; 
	@{$back{'aa2antiC'}{'Val'}} = qw/AAC GAC CAC TAC/; 
	@{$back{'aa2antiC'}{'Ser'}} = qw/AGA GGA CGA TGA ACT GCT/; 
	@{$back{'aa2antiC'}{'Arg'}} = qw/ACG GCG CCG TCG CCT TCT/; 
	@{$back{'aa2antiC'}{'Leu'}} = qw/AAG GAG CAG TAG CAA TAA/; 
	@{$back{'aa2antiC'}{'Phe'}} = qw/AAA GAA/; 
	@{$back{'aa2antiC'}{'Asn'}} = qw/ATT GTT/; 
	@{$back{'aa2antiC'}{'Lys'}} = qw/CTT TTT/; 
	@{$back{'aa2antiC'}{'Asp'}} = qw/ATC GTC/; 
	@{$back{'aa2antiC'}{'Glu'}} = qw/CTC TTC/; 
	@{$back{'aa2antiC'}{'His'}} = qw/ATG GTG/; 
	@{$back{'aa2antiC'}{'Gln'}} = qw/CTG TTG/; 
	@{$back{'aa2antiC'}{'Ile'}} = qw/AAT GAT TAT/; # In stats: 'AAT GAT CAT TAT'; 
	@{$back{'aa2antiC'}{'Met'}} = qw/CAT/; # Also iMet; 
	@{$back{'aa2antiC'}{'Tyr'}} = qw/ATA GTA/; 
	@{$back{'aa2antiC'}{'Supres'}} = qw/CTA TTA/; # Stop codon. In stats : 'CTA TTA TCA'
	@{$back{'aa2antiC'}{'Cys'}} = qw/ACA GCA/; 
	@{$back{'aa2antiC'}{'Trp'}} = qw/CCA/; 
	@{$back{'aa2antiC'}{'SelCys'}} = qw/TCA/; # Stop codon. 
	@{$back{'FourBox'}} = qw/Ala Gly Pro Thr Val/; 
	@{$back{'SixBox'}}  = qw/Ser Arg Leu/; 
	@{$back{'TwoBox'}}  = qw/Phe Asn Lys Asp Glu His Gln/; 
	@{$back{'OtherBox'}} = qw/Ile Met Tyr Supres Cys Trp SelCys/; 
	for my $aa (sort keys %{$back{'aa2antiC'}}) {
		for my $ac (@{$back{'aa2antiC'}{$aa}}) {
			defined $back{'antiC2aa'}{$ac} and do { &tsmsg("[Wrn] Skip repeated anticodon [$ac]\n"); next; }; 
			$back{'antiC2aa'}{$ac} = $aa; 
			$back{'antiC2C'}{$ac} = reverse($ac); 
			$back{'antiC2C'}{$ac} =~ tr/acgturykmbvdhACGTURYKMBVDHwWsSnN/tgcaayrmkvbhdTGCAAYRMKVBHDwWsSnN/; 
			$back{'C2antiC'}{ $back{'antiC2C'}{$ac} } = $ac; 
		}
	}
	return(%back); 
}# iso_codon_table 


sub init_ofiles {
	$ofn{'boxTbl'} = "$opts{'opref'}.box"; 
	&fileSunhh::write2file($ofn{'boxTbl'}, '', '>'); 
	&fileSunhh::write2file($ofn{'boxTbl'}, join("\t", qw/BoxType Isotype total/, map { "antiC_$_" } (1..6))."\n", '>'); 
	if (defined $opts{'inCds'}) {
		$ofn{'codonTbl'} = "$opts{'opref'}.codon"; 
		&fileSunhh::write2file($ofn{'codonTbl'}, join("\t", qw/BoxType Isotype total/, map { "codon_$_" } (1..6))."\n", '>'); 
		$ofn{'codon2gene'} = "$opts{'opref'}.codon2gene"; 
		&fileSunhh::write2file($ofn{'codon2gene'}, join("\t", qw/Isotype anticodon codon genNum genes/)."\n", '>'); 
	}

}#

sub load_cds {
	%cds = %{$fs_obj->save_seq_to_hash('faFile'=>$opts{'inCds'})}; 
	for my $k1 (keys %cds) {
		my $frame = 1; 
		$cds{$k1}{'seq'} =~ s!\s!!g; 
		$cds{$k1}{'len'} = length($cds{$k1}{'seq'}); 
		$cds{$k1}{'head'} =~ m!\[frame=([+-]?\d+)\]!i and $frame = $1; 
		$frame =~ m!^[+-]?(1|2|3)$! or do { &tsmsg("[Wrn] Skip bad frame information [$frame]\n"); $frame = 1; }; 
		$cds{$k1}{'frame'} = $frame; 
		my $t_seq; 
		if ( $frame == 1 ) {
			$t_seq = $cds{$k1}{'seq'}; 
		} elsif ( $frame == 2 ) {
			$t_seq = substr($cds{$k1}{'seq'}, 2); 
		} elsif ( $frame == 3 ) {
			$t_seq = substr($cds{$k1}{'seq'}, 1); 
		} else {
			&stopErr("[Err] Failed to parse bad frame [$frame]\n"); 
		}
		my $t_len = length($t_seq); 
		for (my $i=0; $i<$t_len; $i+=3) {
			$i+3 <= $t_len or last; 
			my $bbb = uc( substr($t_seq, $i, 3) ); 
			$bbb =~ m!^[ATGC]{3}$! or next; 
			$codonInCDS{'codonCnt'}{$bbb} ++; 
			$codonInCDS{'codon2gene'}{$bbb}{$k1}++; 
			$codonInCDS{'codonAllCnt'}++; 
		}
	}
	for my $bbb (keys %{$codonInCDS{'codonCnt'}}) {
		$codonInCDS{'codonRatio'}{$bbb} = $codonInCDS{'codonCnt'}{$bbb} / $codonInCDS{'codonAllCnt'}; 
	}
}# load_cds() 

