#!/usr/bin/perl -w
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 

# [Sunhh@whale 00rawGenom]$ head -3 P1Genom_V1p2.prot.fa P1Genom_V1p2.prot_chr.gff3
# ==> P1Genom_V1p2.prot.fa <==
# >Cma_029951
# MAWINAVLQRPPTCCLLPSLSLLDRKMRIAKAFPPREVQTSSEREKRIETRSQIPSESKV
# KAERCELLVSLIPLSLPFPFASNSSTYLEKLSAYECGFDPSGDARSRFDIRFYLVSILFI
#
# ==> P1Genom_V1p2.prot_chr.gff3 <==
# Cma_Chr00       maker   gene    55660   57727   .       +       .       ID=Cma_031679-gene;Name=Cma_031679-gene
# Cma_Chr00       maker   mRNA    55660   57727   .       +       .       ID=Cma_031679;Parent=Cma_031679-gene;Name=Cma_031679;_AED=0.87;_eAED=0.90;_QI=0|0|0|0.5|0|0|6|0|241
# Cma_Chr00       maker   exon    55660   55677   .       +       .       ID=Cma_031679:exon:21;Parent=Cma_031679

!@ARGV and die "perl $0 P1Genom_V1p2.prot.fa P1Genom_V1p2.prot_chr.gff3 outPref badIDs\n"; 

$ARGV[2] //= 'OUT0'; # This must be exactly four characters length! 
$ARGV[3] //= ''; 

$ARGV[2] = sprintf("%4.4s", $ARGV[2]); 
$ARGV[2] =~ s/\s/_/g; 

my %bad_chrID; 
if ($ARGV[3] ne '') {
	for my $ta (split(/,/, $ARGV[3])) {
		$ta =~ s/\s//g; 
		$bad_chrID{$ta} = 1; 
	}
}

my %prot_fa = %{ $fs_obj->save_seq_to_hash( 'faFile'=>$ARGV[0] ) }; 

my %h; 
$h{'chrN'} = 0; 
$h{'IDf_all'} = 0; 
$h{'IDg_all'} = 0; 
my $fh = &openFH( $ARGV[1], '<' ); 
my $o1 = &openFH( "$ARGV[2].ch",  '>' ); 
my $o2 = &openFH( "$ARGV[2].def", '>' ); 
my $o3 = &openFH( "$ARGV[2].prt", '>' ); 
print {$o2} join( "\t", qw(type    name            chr     start   end     strand  sens    IDg/chr IDg/all IDf/all) )."\n"; 
while (&wantLineC($fh)) {
	my @ta = &splitL("\t", $_); 
	defined $bad_chrID{$ta[0]} and next; 
	if ( $ta[2] eq 'mRNA' ) {
		$ta[8] =~ m/^ID=([^\s;]+)/ or die "[Err] $_\n"; 
		my $mID = $1; 
		defined $prot_fa{ $mID } or do { &tsmsg("[Wrn] Missing seq for mID=$mID\n"); next; }; 
		defined $h{'mID'}{$mID} and next; 
		defined $h{'chr'}{$ta[0]} or do { $h{'chrN'}++; $h{'chr'}{$ta[0]} = sprintf("%05d", $h{'chrN'}); $h{'IDg_chr'}{$ta[0]} = 0; }; 
		$h{'IDf_all'} ++; 
		$h{'IDg_all'} ++; 
		$h{'IDg_chr'}{$ta[0]} ++; 
		$h{'IDf_all_chr'}{$ta[0]} = $h{'IDf_all'}; 
		my $tag_tf = 't'; $ta[6] eq '-' and $tag_tf = 'f'; 
		my @tb = (
		  $h{'chr'}{$ta[0]}, 
		  $ta[3], 
		  $ta[4], 
		  $ta[6], 
		  $tag_tf, 
		  sprintf("%05d", $h{'IDg_chr'}{$ta[0]}), 
		  sprintf("%05d", $h{'IDg_all'}), 
		  sprintf("%05d", $h{'IDf_all'})
		); 
		print {$o2} join( "\t", 'gene', $mID, @tb )."\n"; 
		my $ss = $prot_fa{ $mID }{'seq'}; 
		$ss =~ s/\s//sg; $ss =~ s/\*+$//; $ss =~ s/(.{60})/$1\n/g; chomp($ss); 
		print {$o3} join('', join("\t", ">$prot_fa{$mID}{'key'}", @tb)."\n" , $ss )."\n"; 
	}
}
close($fh); 
close($o2); 
close($o3); 
my @tks = sort { $h{'chr'}{$a} cmp $h{'chr'}{$b} } keys %{$h{'chr'}}; 
for ( @tks ) {
	$h{'chr'}{$_} += 0; 
}
print {$o1} join("\t", @{$h{'chr'}}{@tks})."\n"; 
print {$o1} join("\t", @{$h{'IDf_all_chr'}}{@tks})."\n"; 
print {$o1} join("\t", @{$h{'IDf_all_chr'}}{@tks})."\n"; 

# [Sunhh@whale 01Genomes]$ head ZYRO.def
# type    name            chr     start   end     strand  sens    IDg/chr IDg/all IDf/all
# pseudoG ZYRO0A00110g    001     1150    2192    +       t       00000   00000   00001
#
# gene    ZYRO0A00132g    001     3142    4002    -       f       00001   00001   00002
# gene    ZYRO0A00154g    001     4324    5022    +       t       00002   00002   00003
# pseudoG ZYRO0A00176g    001     6051    6665    +       t       00002   00002   00004
# gene    ZYRO0A00198g    001     7771    8634    -       f       00003   00003   00005
# centromere      centre_A        001     369077  369243  -       f       00185   00185   00202

# >ZYRO0A00132g   001     3142    4002    -       f       00001   00001   00002
# MSIISDESEMDLEKCVQTFEEEEIELPEDSSGGKLTYFLAPYFFEKKWISDIIFPLFGAL
# LVLLIRLAYFSLRDYFDTQLEMKPSDIFSPLWMIIHYLLCFTTCFVNFIYHRKNFNIETK
# ALQKLLKEVAEMDLTGDPAAWQRIASRVNHFSEEGGHHYPLFYSGEHCMRFFVREIVKPI

# [Sunhh@whale 01Genomes]$ more ZYRO.ch
# # ChrName:        A       B       C       D       E       F       G
# # centromere|end: 202     1088    1743    2666    3247    3727    4855
# # End_IDf/all:    634     1400    2240    3051    3531    4377    5399

