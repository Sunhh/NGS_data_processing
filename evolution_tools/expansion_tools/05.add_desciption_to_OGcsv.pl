#!/usr/bin/perl
# ID	arabidopsis	cucumber	grape	maize	maxima	melon	moschata	papaya	populus	rice	tomato	watermelon
# OG0000036	arabidopsis.ATCG00860.1, arabidopsis.ATCG01280.1	cucumber.Csa2M059780.1	grape.GSVIVT01004896001, grape.GSVIVT01011431001, grape.GSVIVT01029119001, grape.GSVIVT01029605001, grape.GSVIVT01029862001	maize.GRMZM2G029465_P01, maize.GRMZM2G065016_P01, maize.GRMZM2G137382_P01, maize.GRMZM2G325417_P01, maize.GRMZM2G359961_P01, maize.GRMZM2G474511_P01, maize.GRMZM2G700821_P01, maize.GRMZM5G803454_P01, maize.GRMZM5G848620_P01, maize.GRMZM5G857367_P01, maize.GRMZM5G882971_P01	maxima.Cma_007911, maxima.Cma_012157, maxima.Cma_017217, maxima.Cma_019844, maxima.Cma_019846, maxima.Cma_026429, maxima.Cma_027423, maxima.Cma_028603, maxima.Cma_029031, maxima.Cma_029621, maxima.Cma_029622, maxima.Cma_029624, maxima.Cma_029828, maxima.Cma_029830, maxima.Cma_029831, maxima.Cma_029832, maxima.Cma_029833, maxima.Cma_029946, maxima.Cma_029947, maxima.Cma_031059, maxima.Cma_031081, maxima.Cma_031088	melon.MELO3C000900P1, melon.MELO3C001174P1, melon.MELO3C027250P1	moschata.Cmo_000196, moschata.Cmo_006745, moschata.Cmo_006747, moschata.Cmo_006748, moschata.Cmo_006749, moschata.Cmo_011501, moschata.Cmo_022638, moschata.Cmo_029666, moschata.Cmo_029667, moschata.Cmo_031078	papaya.evm.TU.contig_42550.3, papaya.evm.model.supercontig_10.247, papaya.evm.model.supercontig_30.80, papaya.evm.model.supercontig_3885.1, papaya.evm.model.supercontig_68.6	populus.Potri.001G194800.1, populus.Potri.003G067500.1, populus.Potri.004G043800.1, populus.Potri.004G128100.1, populus.Potri.006G162000.1, populus.Potri.007G054700.1, populus.Potri.007G054800.1, populus.Potri.009G010000.1, populus.Potri.011G075000.1, populus.Potri.013G079300.1, populus.Potri.013G138500.1, populus.Potri.013G138700.1, populus.Potri.013G143400.1, populus.Potri.019G028400.1, populus.Potri.T016500.1		tomato.Solyc00g005080.1.1, tomato.Solyc00g316530.1.1, tomato.Solyc01g007640.2.1, tomato.Solyc01g017060.1.1, tomato.Solyc01g017110.1.1, tomato.Solyc01g017120.1.1, tomato.Solyc01g017130.1.1, tomato.Solyc01g017140.1.1, tomato.Solyc01g017280.1.1, tomato.Solyc01g017290.1.1, tomato.Solyc01g017410.1.1, tomato.Solyc01g020510.1.1, tomato.Solyc01g056460.1.1, tomato.Solyc01g056870.1.1, tomato.Solyc01g080710.1.1, tomato.Solyc02g010620.1.1, tomato.Solyc02g011920.1.1, tomato.Solyc02g011930.1.1, tomato.Solyc02g011940.1.1, tomato.Solyc02g011950.1.1, tomato.Solyc02g011960.1.1, tomato.Solyc02g011970.1.1, tomato.Solyc02g014780.1.1, tomato.Solyc02g014790.1.1, tomato.Solyc02g014800.1.1, tomato.Solyc02g014810.1.1, tomato.Solyc02g014820.1.1, tomato.Solyc02g078420.1.1, tomato.Solyc03g013480.1.1, tomato.Solyc03g013490.1.1, tomato.Solyc03g020090.1.1, tomato.Solyc03g020100.2.1, tomato.Solyc03g059380.1.1, tomato.Solyc03g078010.1.1, tomato.Solyc03g080130.1.1, tomato.Solyc03g093320.1.1, tomato.Solyc03g120370.1.1, tomato.Solyc04g024540.1.1, tomato.Solyc04g024550.1.1, tomato.Solyc04g024560.1.1, tomato.Solyc04g024580.1.1, tomato.Solyc04g026320.1.1, tomato.Solyc04g039760.1.1, tomato.Solyc04g039770.1.1, tomato.Solyc04g039780.1.1, tomato.Solyc04g050840.1.1, tomato.Solyc04g050850.1.1, tomato.Solyc05g026320.1.1, tomato.Solyc05g041680.1.1, tomato.Solyc05g047670.1.1, tomato.Solyc06g016830.1.1, tomato.Solyc06g017930.1.1, tomato.Solyc06g018000.1.1, tomato.Solyc06g032740.1.1, tomato.Solyc06g048990.1.1, tomato.Solyc07g025240.1.1, tomato.Solyc08g014590.1.1, tomato.Solyc08g028670.1.1, tomato.Solyc08g029400.1.1, tomato.Solyc08g029410.1.1, tomato.Solyc08g045620.1.1, tomato.Solyc08g048230.1.1, tomato.Solyc08g062850.1.1, tomato.Solyc09g015480.1.1, tomato.Solyc09g031800.1.1, tomato.Solyc09g031810.1.1, tomato.Solyc09g031820.1.1, tomato.Solyc09g031830.1.1, tomato.Solyc09g031840.1.1, tomato.Solyc09g031850.1.1, tomato.Solyc09g066380.1.1, tomato.Solyc10g012230.1.1, tomato.Solyc10g018980.1.1, tomato.Solyc10g018990.1.1, tomato.Solyc10g019000.1.1, tomato.Solyc10g019010.1.1, tomato.Solyc10g047440.1.1, tomato.Solyc10g078870.1.1, tomato.Solyc11g011430.1.1, tomato.Solyc11g021070.1.1, tomato.Solyc11g021080.1.1, tomato.Solyc11g030820.1.1, tomato.Solyc11g040080.1.1, tomato.Solyc11g056520.1.1, tomato.Solyc11g062050.1.1, tomato.Solyc11g065050.1.1, tomato.Solyc11g066570.1.1, tomato.Solyc12g019590.1.1, tomato.Solyc12g033010.1.1, tomato.Solyc12g035760.1.1, tomato.Solyc12g035770.1.1, tomato.Solyc12g035780.1.1, tomato.Solyc12g035790.1.1, tomato.Solyc12g035970.1.1, tomato.Solyc12g035980.1.1, tomato.Solyc12g035990.1.1, tomato.Solyc12g036000.1.1, tomato.Solyc12g062470.1.1, tomato.Solyc12g070220.1.1	watermelon.Cla001286, watermelon.Cla013972, watermelon.Cla021817
# OG0000041	arabidopsis.AT1G50060.1, arabidopsis.AT2G14580.1, arabidopsis.AT2G14610.1, arabidopsis.AT2G19990.1, arabidopsis.AT3G09590.1, arabidopsis.AT3G19690.1, arabidopsis.AT4G25790.1, arabidopsis.AT4G30320.1, arabidopsis.AT4G33710.1, arabidopsis.AT4G33720.1, arabidopsis.AT4G33730.1, arabidopsis.AT5G02730.1, arabidopsis.AT5G26130.1, arabidopsis.AT5G57625.1	cucumber.Csa7M070230.1	grape.GSVIVT01029124001, grape.GSVIVT01036993001, grape.GSVIVT01036997001, grape.GSVIVT01037000001, grape.GSVIVT01037005001, grape.GSVIVT01037008001, grape.GSVIVT01037011001, grape.GSVIVT01037015001, grape.GSVIVT01038540001	maize.AC205274.3_FGP001, maize.AC211357.4_FGP002, maize.GRMZM2G053493_P01, maize.GRMZM2G304442_P01, maize.GRMZM2G437187_P01, maize.GRMZM2G456997_P01, maize.GRMZM2G465226_P01, maize.GRMZM5G852886_P01	maxima.Cma_000421, maxima.Cma_011460, maxima.Cma_012165, maxima.Cma_012166, maxima.Cma_012167, maxima.Cma_028658, maxima.Cma_028661, maxima.Cma_028662, maxima.Cma_028663, maxima.Cma_028664, maxima.Cma_028665, maxima.Cma_028666, maxima.Cma_031409, maxima.Cma_031480, maxima.Cma_031481, maxima.Cma_031482, maxima.Cma_031483, maxima.Cma_031486, maxima.Cma_031487, maxima.Cma_031488, maxima.Cma_031654, maxima.Cma_031655, maxima.Cma_031656, maxima.Cma_031657, maxima.Cma_031662, maxima.Cma_031663, maxima.Cma_031664, maxima.Cma_031665	melon.MELO3C017496P1, melon.MELO3C018034P1, melon.MELO3C018536P1, melon.MELO3C018538P1, melon.MELO3C018539P1, melon.MELO3C018540P1, melon.MELO3C018544P1, melon.MELO3C018545P1, melon.MELO3C018546P1, melon.MELO3C018547P1	moschata.Cmo_001609, moschata.Cmo_001610, moschata.Cmo_001611, moschata.Cmo_001612, moschata.Cmo_001613, moschata.Cmo_001614, moschata.Cmo_001615, moschata.Cmo_001616, moschata.Cmo_001617, moschata.Cmo_001618, moschata.Cmo_001619, moschata.Cmo_003052, moschata.Cmo_009312, moschata.Cmo_009313, moschata.Cmo_016336, moschata.Cmo_016337, moschata.Cmo_016338, moschata.Cmo_016341, moschata.Cmo_016342, moschata.Cmo_016343, moschata.Cmo_016344, moschata.Cmo_016345, moschata.Cmo_016346, moschata.Cmo_016348, moschata.Cmo_016350, moschata.Cmo_021622, moschata.Cmo_024676, moschata.Cmo_028545, moschata.Cmo_028546, moschata.Cmo_028547, moschata.Cmo_030017, moschata.Cmo_030018, moschata.Cmo_030019, moschata.Cmo_030021	papaya.evm.TU.contig_27001.2, papaya.evm.model.supercontig_20.79, papaya.evm.model.supercontig_20.80, papaya.evm.model.supercontig_95.71	populus.Potri.001G288400.1, populus.Potri.001G288600.1, populus.Potri.006G171300.1, populus.Potri.009G082800.1, populus.Potri.009G082900.1, populus.Potri.009G083000.1, populus.Potri.009G083100.1, populus.Potri.009G083300.1, populus.Potri.009G083600.1, populus.Potri.T093500.1, populus.Potri.T131400.1, populus.Potri.T131500.1	rice.LOC_Os01g28450.1, rice.LOC_Os01g28500.1, rice.LOC_Os02g27300.1, rice.LOC_Os02g54530.1, rice.LOC_Os02g54540.1, rice.LOC_Os06g24290.1, rice.LOC_Os07g03279.1, rice.LOC_Os07g03288.1, rice.LOC_Os07g03319.1, rice.LOC_Os07g03368.1, rice.LOC_Os07g03377.1, rice.LOC_Os07g03409.1, rice.LOC_Os07g03458.1, rice.LOC_Os07g03467.1, rice.LOC_Os07g03499.1, rice.LOC_Os07g03580.1, rice.LOC_Os07g03590.1, rice.LOC_Os07g03600.1, rice.LOC_Os07g03610.1, rice.LOC_Os07g03620.1, rice.LOC_Os07g03680.1, rice.LOC_Os07g03690.1, rice.LOC_Os07g03710.1, rice.LOC_Os07g03730.1, rice.LOC_Os07g03740.1, rice.LOC_Os07g03750.1, rice.LOC_Os07g14030.1, rice.LOC_Os10g11500.1	tomato.Solyc00g174340.1.1, tomato.Solyc01g106600.2.1, tomato.Solyc01g106610.2.1, tomato.Solyc01g106620.2.1, tomato.Solyc01g106640.2.1, tomato.Solyc07g006710.1.1, tomato.Solyc09g006010.2.1, tomato.Solyc09g007010.1.1, tomato.Solyc09g007020.1.1	watermelon.Cla001246, watermelon.Cla001621, watermelon.Cla001623, watermelon.Cla001624, watermelon.Cla001628, watermelon.Cla001629, watermelon.Cla001630, watermelon.Cla001631, watermelon.Cla006468, watermelon.Cla006469, watermelon.Cla018801
# OG0000042	arabidopsis.AT1G12090.1, arabidopsis.AT1G12100.1, arabidopsis.AT1G62510.1, arabidopsis.AT2G45180.1, arabidopsis.AT4G12470.1, arabidopsis.AT4G12480.1, arabidopsis.AT4G12490.1, arabidopsis.AT4G12500.1, arabidopsis.AT4G12510.1, arabidopsis.AT4G12520.1, arabidopsis.AT4G12545.1, arabidopsis.AT4G12550.1, arabidopsis.AT4G22460.1, arabidopsis.AT5G46890.1, arabidopsis.AT5G46900.1	cucumber.Csa3M021150.1, cucumber.Csa3M027210.1, cucumber.Csa5M161310.1, cucumber.Csa5M161340.1, cucumber.Csa5M161350.1, cucumber.Csa5M161850.1, cucumber.Csa5M161860.1, cucumber.Csa5M161870.1, cucumber.Csa5M161880.1, cucumber.Csa5M161890.1, cucumber.Csa5M161900.1	grape.GSVIVT01001296001, grape.GSVIVT01001297001, grape.GSVIVT01001299001	maize.AC155352.2_FGP010, maize.AC234161.1_FGP001, maize.AC234161.1_FGP002, maize.GRMZM2G027167_P01, maize.GRMZM2G037255_P01, maize.GRMZM2G091534_P01, maize.GRMZM2G094639_P01, maize.GRMZM2G104945_P01, maize.GRMZM2G106324_P01, maize.GRMZM2G162276_P01, maize.GRMZM2G163582_P01, maize.GRMZM2G351505_P01, maize.GRMZM2G379898_P01, maize.GRMZM2G391272_P01, maize.GRMZM2G391286_P01, maize.GRMZM2G398807_P01, maize.GRMZM2G406313_P01, maize.GRMZM2G410338_P01, maize.GRMZM2G429000_P01, maize.GRMZM2G477685_P01	maxima.Cma_008932, maxima.Cma_016732, maxima.Cma_022910, maxima.Cma_022912, maxima.Cma_022914, maxima.Cma_027959, maxima.Cma_027960, maxima.Cma_027961, maxima.Cma_027962, maxima.Cma_027963, maxima.Cma_027965, maxima.Cma_027966, maxima.Cma_027968, maxima.Cma_027969, maxima.Cma_030529, maxima.Cma_030530	melon.MELO3C005540P1, melon.MELO3C005541P1, melon.MELO3C005542P1, melon.MELO3C005543P1, melon.MELO3C005544P1, melon.MELO3C005545P1, melon.MELO3C005546P1, melon.MELO3C005547P1, melon.MELO3C005548P1, melon.MELO3C005549P1, melon.MELO3C005550P1, melon.MELO3C013952P1, melon.MELO3C013954P1, melon.MELO3C013958P1	moschata.Cmo_000702, moschata.Cmo_000703, moschata.Cmo_000704, moschata.Cmo_000705, moschata.Cmo_000706, moschata.Cmo_000707, moschata.Cmo_000814, moschata.Cmo_000816, moschata.Cmo_000817, moschata.Cmo_020029, moschata.Cmo_024844, moschata.Cmo_024845, moschata.Cmo_024846, moschata.Cmo_024847, moschata.Cmo_024848, moschata.Cmo_024849, moschata.Cmo_024850, moschata.Cmo_024851, moschata.Cmo_024853, moschata.Cmo_024854	papaya.evm.model.supercontig_276.2, papaya.evm.model.supercontig_276.3, papaya.evm.model.supercontig_276.4, papaya.evm.model.supercontig_475.2	populus.Potri.001G121800.1, populus.Potri.001G121900.1, populus.Potri.001G158400.1, populus.Potri.003G111400.1, populus.Potri.006G256100.1, populus.Potri.013G128800.1, populus.Potri.018G025900.1, populus.Potri.T145400.1	rice.LOC_Os02g44310.1, rice.LOC_Os02g44320.1, rice.LOC_Os03g01300.1, rice.LOC_Os03g01320.1, rice.LOC_Os03g58670.1, rice.LOC_Os04g46810.1, rice.LOC_Os04g46820.1, rice.LOC_Os04g46830.1, rice.LOC_Os04g52260.1, rice.LOC_Os06g01580.1, rice.LOC_Os08g27678.1, rice.LOC_Os10g09920.1, rice.LOC_Os10g20830.1, rice.LOC_Os10g20840.1, rice.LOC_Os10g20860.1, rice.LOC_Os10g20880.1, rice.LOC_Os10g20890.1, rice.LOC_Os10g40420.1, rice.LOC_Os10g40430.1, rice.LOC_Os10g40440.1, rice.LOC_Os10g40460.1, rice.LOC_Os10g40470.1, rice.LOC_Os10g40480.1, rice.LOC_Os10g40510.1, rice.LOC_Os10g40520.1, rice.LOC_Os10g40530.1, rice.LOC_Os10g40614.1	tomato.Solyc03g083990.1.1, tomato.Solyc03g090990.1.1, tomato.Solyc03g091000.1.1, tomato.Solyc03g091010.1.1, tomato.Solyc03g091020.1.1, tomato.Solyc03g091030.1.1, tomato.Solyc03g091040.1.1, tomato.Solyc03g093040.2.1, tomato.Solyc03g093050.1.1, tomato.Solyc03g093060.1.1, tomato.Solyc03g093070.1.1, tomato.Solyc06g065970.1.1, tomato.Solyc08g005960.1.1, tomato.Solyc08g074480.1.1, tomato.Solyc08g078900.1.1, tomato.Solyc08g078910.1.1, tomato.Solyc08g078930.1.1, tomato.Solyc08g078940.1.1, tomato.Solyc12g014620.1.1	watermelon.Cla005657, watermelon.Cla005659, watermelon.Cla005662, watermelon.Cla011569, watermelon.Cla011570, watermelon.Cla011571, watermelon.Cla011572, watermelon.Cla011573, watermelon.Cla011574
use strict; 
use warnings; 
use Getopt::Long; 
use LogInforSunhh; 
use fileSunhh; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"tab_desc:s@", 
	"noTrim!", 
); 

my $help_txt = <<HH; 

perl $0 -tab_desc P1_ahrd.2c   Orthogroups.csv > Orthogroups.csv.desc

-noTrim     [Bool]

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my %g2desc; 
for my $fn (@{$opts{'tab_desc'}}) {
	for my $tr1 (&fileSunhh::load_tabFile($fn)) {
		push(@{$g2desc{$tr1->[0]}}, $tr1->[1]); 
	}
}
sub arr2txt {
	my ($hr) = @_; 
	my %back; 
	for my $k1 (sort keys %$hr) {
		my %cnt; 
		my $r=0; 
		for my $d1 (@{$hr->{$k1}}) {
			my @a2 = &splitL("_;;;_", $d1); 
			for my $a3 (@a2) {
				$r++; 
				if ($a3 =~ s! \((\d+)\)_$!!) {
					$cnt{$a3}{'n'} += $1; 
					$cnt{$a3}{'r'} //= $r; 
				} else {
					$cnt{$a3}{'n'} ++; 
					$cnt{$a3}{'r'} //= $r; 
				}
			}
		}
		$back{$k1} = join('_;;;_', map { "${_} ($cnt{$_}{'n'})_" } (sort { $cnt{$b}{'n'}<=>$cnt{$a}{'n'} || $cnt{$a}{'r'} <=> $cnt{$b}{'r'} } keys %cnt)); 
	}
	return(\%back); 
}# arr2txt() 
%g2desc = %{ &arr2txt( \%g2desc ) }; 

while (<>) {
	chomp; 
	my @ta = &splitL("\t", $_); 
	if ($ta[0] eq '' or $ta[0] eq 'ID') {
		print join("\t", @ta, 'Description')."\n"; 
		next; 
	}
	my @desc_arr; 
	for my $t1 (@ta[1..$#ta]) {
		for my $t2 (&splitL(",", $t1)) {
			$t2 =~ s!^\s+!!; $t2 =~ s!\s+$!!; 
			$t2 eq '' and next; 
			if (defined $g2desc{$t2}) {
				push(@desc_arr, $g2desc{$t2} ); 
			} elsif (!$opts{'noTrim'}) {
				my $t3 = $t2; 
				$t3 =~ s!^([^\s.]+)\.!!; 
				defined $g2desc{$t3} and push(@desc_arr, $g2desc{$t3}); 
			}
		}
	}
	my $t4 = &arr2txt( { 'xx' => \@desc_arr} ); 
	print join("\t", @ta, $t4->{'xx'})."\n"; 
}


