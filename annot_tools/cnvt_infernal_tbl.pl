#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"outBed!", 
	"outGff!", 
	"noTrunc!", 
); 


my $usage = <<HH; 

perl $0 P1Genom_V1.scf.fa.rRNA.tbl > P1Genom_V1.scf.fa.rRNA.tbl.tab

-outBed
-outGff

-noTrunc

HH

-t and !@ARGV and &LogInforSunhh::usage($usage); 
$opts{'help'} and &LogInforSunhh::usage($usage);

my @lines; 
while (<>) {
	chomp; 
	push(@lines, $_); 
}
if ($lines[1] =~ m!^#\-!) {
	my @hd; 
	while ($lines[1] =~ s!^([#\-]+\s*)!!) {
		my $eL = length($1); 
		push(@hd, substr($lines[0], 0, $eL)); 
		substr($lines[0], 0, $eL) = ''; 
		$hd[-1] =~ s/^\s+|\s+$//g; 
		$hd[-1] =~ s/\s/_/g; 
	}
	$lines[1] eq '' or die "[Err] Bad header line 2:$lines[1]\n"; 
	$lines[0] = [ @hd ]; 
	splice(@lines, 1,1); 
} else {
	my @hd; 
	for my $k1 (qw/target_name query_name mdl_from mdl_to seq_from seq_to description_of_target/) {
		my $k2 = $k1; 
		$k2 =~ s!_! !g; 
		$lines[1] =~ s!$k2!$k1!g; 
	}
	$lines[0] = [ split(/\s+/, $lines[0]) ]; 
}
$lines[0][0] =~ s/^#//; 

#target_name	accession	query_name	accession	mdl	mdl_from	mdl_to	seq_from	seq_to	strand	trunc	pass	gc	bias	score	E-value	inc	description_of_target

##target name         accession query name           accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of
##------------------- --------- -------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- ---------------
#Cmo_Scf00021         -         5S_rRNA              RF00001    cm        1      119  4452908  4453026      +    no    1 0.55   0.0  126.3     7e-27 !   -

if ( $opts{'outBed'} ) {
	&out_file(\@lines, 'outBed'); 
} elsif ( $opts{'outGff'} ) {
	&out_file(\@lines, 'outGff'); 
} else {
	&out_file(\@lines, 'outTbl'); 
}

sub out_file {
	my ($lines_aref, $outType) = @_; 
	$outType //= 'outTbl'; 
	
	my $hd = shift(@$lines_aref); 
	if ($outType eq 'outTbl') {
		print STDOUT join("\t", @$hd)."\n"; 
	}

	my %h; 
	for my $tl (@{$lines_aref}) {
		$tl =~ m/^\s*(#|$)/ and next; 
		my @ta = split(/\s+/, $tl); 
		$ta[7] > $ta[8] and @ta[7,8] = @ta[8,7]; 
		$opts{'noTrunc'} and do { $ta[10] =~ m/^no$/i or next; }; 
		if ($outType eq 'outTbl') {
			print STDOUT join("\t", @ta)."\n"; 
		} elsif ($outType eq 'outBed') {
			print STDOUT join("\t", $ta[0], $ta[7]-1, $ta[8], '.', 0, $ta[9])."\n"; 
		} elsif ($outType eq 'outGff') {
			my $s1 = 'p'; $ta[9] eq '-' and $s1 = 'm'; 
			my $id = "$ta[0]__$ta[7]__$ta[8]__$s1"; 
			my $new_id = $id; 
			my $cnt = 0; 
			while (defined $h{$new_id}) {
				$cnt++; 
				$new_id = "${id}_n$cnt"; 
				$cnt >= 1e6 and &stopErr("[Err] Something maybe wrong for $new_id. \n"); 
			}
			$id = $new_id; 
			print STDOUT join("\t", $ta[0], 'infernal', 'rRNA', $ta[7], $ta[8], $ta[14], $ta[9], '.', "ID=$id; bias=$ta[13]; eval=$ta[15]; inc=$ta[16]; pass=$ta[11]; trunc=\"$ta[10]\"; query=\"$ta[2]\"")."\n"; 
		}
	}
}




