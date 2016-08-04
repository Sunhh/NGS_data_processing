#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"scf_gff:s", 
	"agp:s", 
	"noSort!", 
	"help!", 
); 

$opts{'help'} and die "perl $0 -scf_gff msc.gff -agp scf2chr.agp [-noSort]\n"; 

( defined $opts{'scf_gff'} and defined $opts{'agp'} ) or die "perl $0 -h\n"; 

my (%ctg2scf) = %{ &fileSunhh::load_agpFile( $opts{'agp'} ) }; 
my %scf_info; 
my $fh1 = &openFH( $opts{'scf_gff'}, '<' ); 
while (<$fh1>) {
	chomp; 
	my @ta = &splitL("\t", $_); 
	push(@{$scf_info{'line'}{$ta[0]}}, [$. , @ta]); 
	$scf_info{'order'}{$ta[0]} //= $.; 
}
close($fh1); 

if ( $opts{'noSort'} ) {
	for my $a1 (sort { $a->[0] <=> $b->[0] } map { @{$scf_info{'line'}{$_}} } keys %{$scf_info{'line'}}) {
		my ( $b1 ) = $ms_obj->switch_position( 'qry2ref'=>\%ctg2scf , 'qryID' => $a1->[1], 'qryPos' => $a1->[3], 'qryStr'=>'+' ); 
		my ( $b2 ) = $ms_obj->switch_position( 'qry2ref'=>\%ctg2scf , 'qryID' => $a1->[1], 'qryPos' => $a1->[4], 'qryStr'=>'+' ); 
		( defined $b1 and defined $b2 ) or die "failed to locate [@$a1]\n"; 
		$a1->[1] = $b1->[0]; 
		$a1->[3] = $b1->[1]; 
		$a1->[4] = $b2->[1]; 
		$a1->[3] > $a1->[4] and ($a1->[3], $a1->[4]) = ($a1->[4], $a1->[3]); 
		print STDOUT join("\t", @{$a1}[ 1 .. $#$a1 ])."\n"; 
	}
} else {
	my %scfP; 
	for my $k1 (keys %{$scf_info{'line'}}) {
		my $a1 = $scf_info{'line'}{$k1}[0]; 
		my ( $b1 ) = $ms_obj->switch_position( 'qry2ref'=>\%ctg2scf , 'qryID' => $a1->[1], 'qryPos' => $a1->[3], 'qryStr'=>'+' ); 
		# my ( $b2 ) = $ms_obj->switch_position( 'qry2ref'=>\%ctg2scf , 'qryID' => $a1->[1], 'qryPos' => $a1->[3], 'qryStr'=>'+' ); 
		# ( defined $b1 and defined $b2 ) or die "failed to locate [@$a1]\n"; 
		defined $b1 or die "failed to locate [@$a1]\n"; 
		$scfP{$k1} = [$b1->[0], $b1->[1], $b1->[2]]; 
		$b1->[2] eq '-' and @{$scf_info{'line'}{$k1}} = reverse(@{$scf_info{'line'}{$k1}}); 
	}
	for my $a1 (map { @{$scf_info{'line'}{$_}} } sort { $scfP{$a}[0] cmp $scfP{$b}[0] || $scfP{$a}[1] <=> $scfP{$b}[1] } keys %{$scf_info{'line'}}) {
		my ( $b1 ) = $ms_obj->switch_position( 'qry2ref'=>\%ctg2scf , 'qryID' => $a1->[1], 'qryPos' => $a1->[3], 'qryStr'=>'+' ); 
		my ( $b2 ) = $ms_obj->switch_position( 'qry2ref'=>\%ctg2scf , 'qryID' => $a1->[1], 'qryPos' => $a1->[4], 'qryStr'=>'+' ); 
		( defined $b1 and defined $b2 ) or die "failed to locate [@$a1]\n"; 
		$a1->[1] = $b1->[0]; 
		$a1->[3] = $b1->[1]; 
		$a1->[4] = $b2->[1]; 
		$a1->[3] > $a1->[4] and ($a1->[3], $a1->[4]) = ($a1->[4], $a1->[3]); 
		print STDOUT join("\t", @{$a1}[ 1 .. $#$a1 ])."\n"; 
	}
}

