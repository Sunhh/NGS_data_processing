#!/usr/bin/perl -w
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"vcf_tab:s", # Input file coming from "vcf_to_tab"; The first three columns are #CHROM , POS and REF . 
	"P1_colN:i", # 3
	"P2_colN:i", # 4
); 

$opts{'P1_colN'} //= 3; 
$opts{'P2_colN'} //= 4; 

my $help_txt = <<HH; 

perl $0 -vcf_tab Itay_Galap_132offs_2parents_geno.tbl 

-help

-P1_colN       [$opts{'P1_colN'}]
-P2_colN       [$opts{'P2_colN'}]

HH

defined $opts{'vcf_tab'} or &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my %glob; 

$glob{'fh_inTab'} = &openFH($opts{'vcf_tab'}, '<'); 
$glob{'need_header'} = 1; 

while (readline($glob{'fh_inTab'})) {
	chomp($_); 
	my @ta = &splitL("\t", $_); 
	if ($ta[0] =~ m!^(#CHROM|chr|chromosome|chrID|CHROM)$!i) {
		# This is header. 
		$glob{'tab_header'} //= [ @ta ]; 
		next; 
	}
	unless ( defined $glob{'offs_colN'} ) {
		for (my $i=3; $i<@ta; $i++) {
			$i == $opts{'P1_colN'} and next; 
			$i == $opts{'P2_colN'} and next; 
			push(@{$glob{'offs_colN'}}, $i); 
		}
	}
	# Get P1 allele : 
	my @p1_al = &tab_allele( $ta[$opts{'P1_colN'}] ); 
	my @p2_al = &tab_allele( $ta[$opts{'P2_colN'}] ); 
	my %pp_al_class = %{ &class_PP_al( \@p1_al, \@p2_al ) }; 
	
	# Class offsprings
	my %cnt; 
	my %genoCnt; 
	for my $i (@{$glob{'offs_colN'}}) {
		$cnt{'total_N'} ++; 
		my @of_al = &tab_allele( $ta[$i] ); 
		if ($#of_al == 0) {
			$genoCnt{"$of_al[0][0]/$of_al[0][0]"} ++; 
		} else {
			$genoCnt{"$of_al[0][0]/$of_al[1][0]"} ++; 
		}
		my $of_class = &class_off_al( \%pp_al_class, \@of_al ); 
		if ( $of_class eq 'miss' ) {
			$cnt{'miss_N'} ++; 
		} elsif ( $of_class eq 'homo_P1_parent' ) {
			$cnt{'homo_P1'} ++; 
		} elsif ( $of_class eq 'homo_P2_parent' ) {
			$cnt{'homo_P2'} ++; 
		} elsif ( $of_class eq 'hete_both_parent' ) {
			$cnt{'hete_PP'} ++; 
		} elsif ( $of_class =~ m!_non_parent$! ) {
			$cnt{'non_PP'} ++; 
		} elsif ( $of_class =~ m!_bad_parent$! ) {
			$cnt{'bad_PP'} ++; 
		} else {
			$cnt{'other'}{ $of_class } ++; 
		}
	}
	$cnt{'other'}{'genoCnt'} = join(':', map { "|$_|=$genoCnt{$_}" } sort keys %genoCnt); 
	if ( $glob{'need_header'} ) {
		$glob{'tab_header'} //= [ 'chr', 'pos', 'ref' ]; 
		$glob{'tab_header'}[$opts{'P1_colN'}] //= 'P1'; $glob{'tab_header'}[$opts{'P1_colN'}] eq '' and $glob{'tab_header'}[$opts{'P1_colN'}] = 'P1'; 
		$glob{'tab_header'}[$opts{'P2_colN'}] //= 'P2'; $glob{'tab_header'}[$opts{'P2_colN'}] eq '' and $glob{'tab_header'}[$opts{'P2_colN'}] = 'P2'; 
		print STDOUT join("\t", @{$glob{'tab_header'}}[0,1,2, $opts{'P1_colN'}, $opts{'P2_colN'}], qw/total miss homo_P1 homo_P2 hete_PP non_PP bad_PP others/)."\n"; 
		$glob{'need_header'} = 0; 
	}
	for my $tk (qw/total_N miss_N homo_P1 homo_P2 hete_PP non_PP bad_PP/) {
		$cnt{$tk} //= 0; 
	}
	$cnt{'other'} //= { 'NA'=>'NA' }; 
	print STDOUT join("\t", 
	  @ta[0,1,2,$opts{'P1_colN'},$opts{'P2_colN'}], 
	  @cnt{qw/total_N miss_N homo_P1 homo_P2 hete_PP non_PP bad_PP/}, 
	  join(";;", map { "$_=$cnt{'other'}{$_}" } sort keys %{$cnt{'other'}})
	)."\n"; 
}
close($glob{'fh_inTab'}); 

=head2 class_off_al( \%parent_al_class, \@off_al )

Return : '(homo|hete)_(non|bad|any|P1|P2|both)_parent' | 'miss' 

  miss                : offspring's genotype is missing. 
  homo_(P1|P2)_parent : homo in offspring, and allele exists in one parent only. 
  hete_(P1|P2)_parent : hete in offspring, and both alleles exist in one parent only. 
  hete_both_parent    : hete in offspring, and one allele comes from one parent, the other one from another parent. 
  ????_non_parent     : Offspring has an allele not existing in any parent. 
  ????_bad_parent     : One of the parent's genotype is missing. 
  ????_any_parent     : At least one allele in offspring is shared by both parents. 

=cut
sub class_off_al {
	my ($pp_alH, $off_alA) = @_; 
	my $class_name = ''; 
	$off_alA->[0][0] eq '.' and return('miss'); 
	for my $t_al (keys %$pp_alH) {
		if ( $pp_alH->{$t_al}{'from'} eq 'UU' ) {
			if ( $#{$off_alA} == 0 ) {
				return('homo_bad_parent'); 
			} elsif ( $#{$off_alA} == 1 ) {
				return('hete_bad_parent')
			} else {
				&stopErr(join('', "[Err] Strange allele : ", join(":", map { @$_ } @$off_alA))."\n"); 
			}
		}
	}
	if ( $#{$off_alA} == 0 ) {
		if ( !(defined $pp_alH->{$off_alA->[0][0]}) ) {
			return('homo_non_parent'); 
		} elsif ( $pp_alH->{ $off_alA->[0][0] }{'from'} eq 'UU' ) {
			return('homo_bad_parent'); 
		} elsif ( $pp_alH->{ $off_alA->[0][0] }{'from'} eq 'PP' ) {
			return('homo_any_parent'); 
		} elsif ( $pp_alH->{ $off_alA->[0][0] }{'from'} eq 'P1' ) {
			return('homo_P1_parent'); 
		} elsif ( $pp_alH->{ $off_alA->[0][0] }{'from'} eq 'P2' ) {
			return('homo_P2_parent'); 
		} else {
			&stopErr("[Err] Unknown allele [$off_alA->[0][0]] with $pp_alH->{$off_alA->[0][0]}{'from'}\n"); 
		}
	} elsif ( $#{$off_alA} == 1 ) {
		defined $pp_alH->{$off_alA->[0][0]} or return('hete_non_parent'); 
		defined $pp_alH->{$off_alA->[1][0]} or return('hete_non_parent'); 
		my $t1 = $pp_alH->{$off_alA->[0][0]}{'from'}; 
		my $t2 = $pp_alH->{$off_alA->[1][0]}{'from'}; 
		if ( $t1 eq $t2 ) {
			$t1 eq 'P1' and return('hete_P1_parent'); 
			$t1 eq 'P2' and return('hete_P2_parent'); 
			$t1 eq 'PP' and return('hete_any_parent'); 
			$t1 eq 'UU' and return('hete_bad_parent'); 
		} elsif ( ( $t1 eq 'P1' and $t2 eq 'P2') or ($t1 eq 'P2' and $t2 eq 'P1') ) {
			return('hete_both_parent'); 
		} elsif ( $t1 eq 'UU' or $t2 eq 'UU' ) {
			return('hete_bad_parent'); 
		} elsif ( $t1 eq 'PP' or $t2 eq 'PP' ) {
			return('hete_any_parent'); 
		} else {
			&stopErr(join('', "[Err] Strange allele : ", join(":", map { @$_ } @$off_alA))."\n"); 
		}
	} else {
		&stopErr(join('', "[Err] Strange allele : ", join(":", map { @$_ } @$off_alA))."\n"); 
	}
	return; 
}# class_off_al() 


=header2 class_PP_al( \@P1_allele , @P2_allele )

Return    : ( \%al_class )

$al_class{ $allele } = { 'from' => 'P1|P2|PP|UU' , 'cnt' => $num_of_ParentAl, 'type' => '(homo|hete)(P1|P2)' or '(homo|hete)P1_(homo|hete)P2' }; 

=cut
sub class_PP_al {
	my ($p1_alA, $p2_alA) = @_; 
	my %back; # '$allele' 
	if ( $p1_alA->[0][0] eq '.' or $p2_alA->[0][0] eq '.' ) {
		for my $ar1 ( @$p1_alA ) {
			$ar1->[0] eq '.' and next; 
			$back{$ar1->[0]}{'from'} = 'UU'; 
			$back{$ar1->[0]}{'cnt'} += $ar1->[1]; 
			$back{$ar1->[0]}{'type'} = ( $ar1->[1] == 2 ) ? 'homoP1' : 'heteP1' ; 
		}
		for my $ar2 ( @$p2_alA ) {
			$ar2->[0] eq '.' and next; 
			$back{$ar2->[0]}{'from'} = 'UU'; 
			$back{$ar2->[0]}{'cnt'} += $ar2->[1]; 
			if ( defined $back{$ar2->[0]} ) {
				$back{$ar2->[0]}{'type'} .= ( ($ar2->[1] == 2) ? '_homoP2' : '_heteP2' ); 
			} else {
				$back{$ar2->[0]}{'type'}  = ( $ar2->[1] == 2 ) ? 'hhhhP1_homoP2' : 'hhhhP1_heteP2'; 
			}
		}
		return(\%back); 
	} 

	for my $ar1 (@$p1_alA) {
		$back{$ar1->[0]}{'from'} = 'P1'; 
		$back{$ar1->[0]}{'cnt'}  = $ar1->[1]; 
		$back{$ar1->[0]}{'type'} = ( $ar1->[1] == 2 ) ? 'homoP1' : 'heteP1' ; 
	}
	for my $ar2 (@$p2_alA) {
		if (defined $back{$ar2->[0]}) {
			$back{$ar2->[0]}{'from'} = 'PP' ; # This means the allele exists in both parents. 
			$back{$ar2->[0]}{'cnt'} += $ar2->[1]; 
			$back{$ar2->[0]}{'type'} .= ( ( $ar2->[1] == 2 ) ? '_homoP2' : '_heteP2' ) ; 
		} else {
			$back{$ar2->[0]}{'from'} = 'P2'; 
			$back{$ar2->[0]}{'cnt'}  = $ar2->[1]; 
			$back{$ar2->[0]}{'type'} = ( $ar2->[1] == 2 ) ? 'homoP2' : 'heteP2' ; 
		}
	}
	return(\%back); 
}# class_PP_al() 

sub tab_allele {
	my $t1 = shift; 
	if ( $t1 eq './.' ) {
		return( ['.', 2] ); 
	} elsif ( $t1 =~ m!^([ATGC\*]+)/([ATGC\*]+)$!i ) {
		my ($a1, $a2) = sort ( uc($1), uc($2) ); 
		if ( $a1 eq $a2 ) {
			return( [$a1, 2] ); 
		} else {
			return( [$a1, 1], [$a2, 1] ); 
		}
	} else {
		&stopErr("[Err] Failed to parse allele [$t1]\n"); 
	}
	&stopErr("[Err] Why here?\n"); 
	return(); 
}# tab_allele ()


