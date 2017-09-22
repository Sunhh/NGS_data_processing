#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

# InFile  Total_size      Total_Rd_num    Mean_Rd_size    Range_Rd_size   PhredCut        Time
# A04_CGATGT_R1.paired.hq 911022366       6454547 141.144276430244        40-151  Phred33 Fri Aug  1 07:50:15 2014
# A04_CGATGT_R2.paired.hq 900417013       6454547 139.501193964503        40-151  Phred33 Fri Aug  1 07:50:24 2014
# A04_CGATGT_R1.single.hq 37172499        383510  96.927065787072 40-151  Phred33 Fri Aug  1 07:37:25 2014
# A04_CGATGT_R2.single.hq 38147318        347569  109.754661664303        40-151  Phred33 Fri Aug  1 07:37:45 2014
# A04_CGATGT_R1.paired    901384698       6389723 141.067883224359        40-151  Phred33 Fri Aug  1 07:49:47 2014
# A04_CGATGT_R2.paired    890837666       6389723 139.417258932821        40-151  Phred33 Fri Aug  1 07:50:10 2014
# A04_CGATGT_R1.single    42437278        420364  100.953644936293        40-151  Phred33 Fri Aug  1 07:37:22 2014
# A04_CGATGT_R2.single    41443368        371126  111.669265963581        40-151  Phred33 Fri Aug  1 07:37:47 2014


my @endSuff = qw/paired.hq single.hq paired single/; 
my %useSuff; 
for (@endSuff) {
	$useSuff{$_} = 1; 
}
print STDOUT join("\t", qw/Prefix HQ_PairedNum	HQ_PairedBp	HQ_SingleNum	HQ_SingleBp	CleanQ20_PairedNum	CleanQ20_PairedBp	CleanQ20_SingleNum	CleanQ20_SingleBp/)."\n"; 

my %infor; 
my @IDs; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[0] eq 'InFile' and next; 
	$ta[0] =~ m!^(\S+)_(R[12])\.(\S+)! or die "$_\n"; 
	my ($curID, $curRN, $curSuff) = ($1, $2, $3); 
	defined $useSuff{$curSuff} or do { &tsmsg("[Err] Unknown suffix [$curSuff]\n"); next; }; 
	defined $infor{$curID} or push( @IDs, $curID ); 
	$infor{$curID}{$curSuff}{$curRN}{rdBp} = $ta[1]; 
	$infor{$curID}{$curSuff}{$curRN}{rdNo} = $ta[2]; 
}

for my $k1 ( @IDs ) {
	my @oo; 
	for my $curSuff ( @endSuff ) {
		defined $infor{$k1}{$curSuff}{R1}{rdNo} or $infor{$k1}{$curSuff}{R1}{rdNo} = 'NA'; 
		defined $infor{$k1}{$curSuff}{R2}{rdNo} or $infor{$k1}{$curSuff}{R2}{rdNo} = 'NA'; 
		defined $infor{$k1}{$curSuff}{R1}{rdBp} or $infor{$k1}{$curSuff}{R1}{rdBp} = 'NA'; 
		defined $infor{$k1}{$curSuff}{R2}{rdBp} or $infor{$k1}{$curSuff}{R2}{rdBp} = 'NA'; 
		push(@oo, "=$infor{$k1}{$curSuff}{R1}{rdNo}+$infor{$k1}{$curSuff}{R2}{rdNo}"); 
		push(@oo, "=$infor{$k1}{$curSuff}{R1}{rdBp}+$infor{$k1}{$curSuff}{R2}{rdBp}"); 
	}
	print STDOUT join("\t", $k1, @oo)."\n"; 
}


