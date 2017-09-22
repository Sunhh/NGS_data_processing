#!/usr/bin/perl
# 2016-12-02 Count common information from .vcf.tab files. 
use strict; 
use warnings; 
use mathSunhh; 
use LogInforSunhh; 
use SNP_tbl; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"cpuN:i",       # 1
	"colN_start:i", # 3
	"noHeader!", 
	"ret_line:s", 
	"rules:s@", 
); 
$opts{'cpuN'} //= 1; 
$opts{'cpuN'} = int($opts{'cpuN'}); 
$opts{'cpuN'} < 1 and $opts{'cpuN'} = 1; 
$opts{'colN_start'} //= 3; 

my @out_type = qw/cnt_indvN cnt_lmiss cnt_alleleTypeN alleleTypeCnt cnt_homo cnt_hete MinorAF hete2N_cnt_lmiss hete2N_cnt_alleleTypeN hete2N_alleleTypeCnt hete2N_cnt_homo hete2N_MinorAF/; 
my @hav_type = (@out_type, qw/MinorAN MajorAN hete2N_MinorAN hete2N_MajorAN/); 



my $help_txt = '#' x 80 . "\n"; 
$help_txt .= <<HH; 

perl $0 in.vcf.tab > in.vcf.tab.l_cnt

 -help
 -cpuN           [$opts{'cpuN'}] 

 -colN_start     [$opts{'colN_start'}]

 -noHeader       [Boolean] If given, the 1st line will be modified too. 

 -ret_line       [filename] Output lines fitting -rules
   -rules         [] Could be multiple times. "type >|<|=|>=|<= value"
                  Example : "cnt_lmiss >= 295"
                  Types that can be used are listed below. 

Format of vcf.tab : 
  chr \\t pos \\t base(ref) \\t sample1    \\t sample2 ...
  c1  \\t 100 \\t A         \\t A/T        \\t G/G
  c2  \\t 10  \\t G         \\t ./.        \\t */*
  c3  \\t 5   \\t T         \\t AATG/AATG  \\t AATG/*
  c4  \\t 8   \\t N         \\t A/*        \\t A/G
  c4  \\t 10  \\t A         \\t ./.        \\t ./.

HH
$help_txt .= "\nTypes that can be used in -rules :\n"; 
for (my $i=0; $i<@hav_type; $i+= 5) {
	my $ei = $i+4; 
	$ei >= $#hav_type and $ei = $#hav_type; 
	$help_txt .= join("\t", '', @hav_type[$i .. $ei])."\n"; 
}
$help_txt .= "#" x 80 . "\n"; 


$opts{'help'} and &LogInforSunhh::usage($help_txt); 
-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

my @InFp = (); 
if ( !@ARGV ) {
	-t or @InFp = (\*STDIN); 
} else {
	for (@ARGV) {
		push( @InFp, &openFH($_, '<') ); 
	}
}
my $pm; 
$opts{'cpuN'} > 1 and $pm = &LogInforSunhh::get_pm( $opts{'cpuN'} ); 

my @rules; 
my $ofh_ret; 
if ( defined $opts{'ret_line'} ) {
	defined $opts{'rules'} or &stopErr("[Err] Need -rules for -ret_line \n"); 
	$ofh_ret = &openFH( $opts{'ret_line'}, '>'); 
	for (@{$opts{'rules'}}) {
		s!^\s+|\s+$!!g;
		m!^(\S+)\s+(\S+)\s+(\S+)$! or die "rules: $_\n"; 
		push(@rules, [$1, $2, $3]); 
	}
}


for my $fh ( @InFp ) {
	my $header_txt = ''; 
	$opts{'noHeader'} or $header_txt = <$fh>; 
	my $header_txt_raw = $header_txt; 
	if ( $opts{'noHeader'} ) {
		my @th = (qw/chr pos/); 
		if ( $opts{'colN_start'} <= 0 ) {
			$header_txt = join("\t", @out_type)."\n"; 
		} elsif ( $opts{'colN_start'} == 1 ) {
			$header_txt = join("\t", $th[0], @out_type)."\n"; 
		} else {
			for ( my $i=0; $i<$opts{'colN_start'}; $i++ ) {
				$th[$i] //= "H$i"; 
			}
			$header_txt = join("\t", @th, @out_type)."\n"; 
		}
		print STDOUT $header_txt; 
	} else {
		chomp($header_txt); 
		my @ta = &splitL( "\t", $header_txt ); 
		my @th = @ta[ 0 .. ($opts{'colN_start'} - 1) ]; 
		$header_txt = join("\t", @th, @out_type)."\n"; 
		print STDOUT $header_txt; 
	} 
	if ( defined $opts{'ret_line'} ) {
		print {$ofh_ret} $header_txt_raw; 
	}

	if ( defined $pm ) {
		$opts{'cpuN'} > 1 or &stopErr("[Err] cpuN not > 1\n"); 
		# separate files
		my $wrk_dir = &fileSunhh::new_tmp_dir( 'create' => 1 ); 
		my @sub_fn = &fileSunhh::dvd_file( $fh, $opts{'cpuN'}, 'keep_order' => 1, 'with_header' => 0, 'sub_pref' => "$wrk_dir/sub_", 'tmpFile' => "$wrk_dir/base_0" ); 
		# process sub-files
		for my $sfn ( @sub_fn ) {
			my $pid = $pm->start and next; 
			open F,'<',"$sfn" or die; 
			open O,'>',"$sfn.o" or die; 
			defined $opts{'ret_line'} and do { open OR,'>',"$sfn.or" or die; }; 
			while (<F>) {
				chomp; 
				my @ta = &splitL("\t", $_); 
				
				my @th = @ta[0 .. ($opts{'colN_start'} - 1)]; 
				my @tb = @ta[ $opts{'colN_start'} .. $#ta ]; 
				my %h = &cnt_type( \@tb ); 
				print O join( "\t", @th, @h{@out_type} )."\n"; 
				if ( defined $opts{'ret_line'} ) {
					my $is_ok = 1; 
					for my $tr1 ( @rules ) {
						&test_rule( \%h, $tr1->[0], $tr1->[1], $tr1->[2] ) or do { $is_ok = 0; last; }; 
					}
					$is_ok == 1 and print OR "$_\n"; 
				}
			}
			close O; 
			defined $opts{'ret_line'} and close OR; 
			close F; 
			$pm->finish; 
		}
		$pm->wait_all_children; 
		# merge sub-files 
		for my $sfn ( @sub_fn ) {
			open F,'<',"$sfn.o" or die; 
			while (<F>) {
				print STDOUT $_; 
			}
			close F; 
			if ( defined $opts{'ret_line'} ) {
				open F,'<',"$sfn.or" or die; 
				while (<F>) {
					print {$ofh_ret} $_; 
				}
				close F; 
			}
		}
		# delete sub-files
		&fileSunhh::_rmtree($wrk_dir); 
	} else {
		# Normal count. 
		while (<$fh>) {
			chomp; 
			my @ta = &splitL( "\t", $_ ); 
			my @th = @ta[0 .. ($opts{'colN_start'} - 1)]; 
			my @tb = @ta[ $opts{'colN_start'} .. $#ta ]; 
			my %h = &cnt_type( \@tb ); 
			print STDOUT join( "\t", @th, @h{@out_type} )."\n"; 
			if ( defined $opts{'ret_line'} ) {
				my $is_ok = 1; 
				for my $tr1 ( @rules ) {
					&test_rule( \%h, $tr1->[0], $tr1->[1], $tr1->[2] ) or do { $is_ok = 0; last; }; 
				}
				$is_ok == 1 and print {$ofh_ret} "$_\n"; 
			}
		}
	}
}

sub cnt_type {
	my ($ar1) = @_; 
	my %cnt_al; 
	my %cnt_geno; 
	my %cnt_indv; 
	for my $b (@$ar1) {
		$cnt_indv{'total'} ++; 
		if ( $b eq './.' ) {
			$cnt_al{'miss'} += 2; 
			$cnt_indv{'miss'} ++; 
			$cnt_geno{'miss'} ++; 
		} elsif ( $b =~ m!^([ATGC\*N]+)/([ATGC\*N]+)$! ) {
			my ($b1, $b2) = ($1, $2); 
			$cnt_al{'any_al'}{$b1} ++; 
			$cnt_al{'any_al'}{$b2} ++; 
			if ( $b1 eq $b2 ) {
				$cnt_geno{'homo'} ++; 
				$cnt_al{'homo_al'}{$b1} += 2; 
			} else {
				$cnt_indv{'hete'} ++; 
				$cnt_geno{'hete'} ++; 
				$cnt_al{'hete_al'}{$b1} ++; 
				$cnt_al{'hete_al'}{$b2} ++; 
			}
		} else {
			&tsmsg("[Err] Unknown format of vcfTab genotype [$b]\n"); 
			# &stopErr("[Err] Stop at @$ar1\n"); 
			$cnt_al{'miss'} += 2; 
			$cnt_indv{'miss'} ++; 
			$cnt_geno{'miss'} ++; 
		}
	}

	for my $tk (qw/homo hete miss/) {
		$cnt_geno{$tk} //= 0; 
	}

	$cnt_al{'miss'} //= 0; 
	for my $tk (qw/any_al homo_al hete_al/) {
		$cnt_al{$tk} //= {}; 
	}

	for my $tk (qw/total miss hete/) {
		$cnt_indv{$tk} //= 0; 
	}

	my %back; 

	$back{'cnt_indvN'}       = $cnt_indv{'total'}; 

	$back{'cnt_lmiss'}       = $cnt_indv{'miss'}; 

	my @srt_allele  = sort { $cnt_al{'any_al'}{$b} <=> $cnt_al{'any_al'}{$a} || $a cmp $b }  keys %{$cnt_al{'any_al'}}; 
	$back{'cnt_alleleTypeN'} = scalar(@srt_allele); 
	if ( $back{'cnt_alleleTypeN'} > 0 ) {
		$back{'alleleTypeCnt'}   = join(";;", map { "${_}_$cnt_al{'any_al'}{$_}" } @srt_allele); 
	} else {
		$back{'alleleTypeCnt'}   = "N_0"; 
	}
	$back{'cnt_homo'}        = $cnt_geno{'homo'}; 
	$back{'cnt_hete'}        = $cnt_geno{'hete'}; 
	my $al_total = &mathSunhh::_sum( map { $cnt_al{'any_al'}{$_} } @srt_allele ); 
	$al_total //= 0; 
	my $al_minor = ( $back{'cnt_alleleTypeN'} > 1 ) ? $cnt_al{'any_al'}{ $srt_allele[1] } : 0 ; 
	# $back{'MinorAF'}         = ( $al_total > 0 ) ? sprintf("%.2f", 100 * $al_minor / $al_total)  : -1 ; 
	$back{'MinorAF'}         = ( $al_total > 0 ) ? (100 * $al_minor / $al_total)  : -1 ; 
	$back{'MinorAN'}         = $al_minor; 
	$back{'MajorAN'}         = ( $back{'cnt_alleleTypeN'} > 0 ) ? $cnt_al{'any_al'}{ $srt_allele[0] } : 0 ; 

	$back{'hete2N_cnt_lmiss'}         = $cnt_indv{'miss'} + $cnt_indv{'hete'}; 
	my @srt_al_homo = sort { $cnt_al{'homo_al'}{$b} <=> $cnt_al{'homo_al'}{$a} || $a cmp $b } keys %{$cnt_al{'homo_al'}}; 
	$back{'hete2N_cnt_alleleTypeN'}   = scalar(@srt_al_homo); 
	if ( $back{'hete2N_cnt_alleleTypeN'} > 0 ) {
		$back{'hete2N_alleleTypeCnt'}     = join(";;", map { "${_}_$cnt_al{'homo_al'}{$_}" } @srt_al_homo); 
	} else {
		$back{'hete2N_alleleTypeCnt'}     = 'N_0'; 
	}
	$back{'hete2N_cnt_homo'}          = $cnt_geno{'homo'}; 
	my $al_total_homo = &mathSunhh::_sum( map { $cnt_al{'homo_al'}{$_} } @srt_al_homo ); 
	$al_total_homo //= 0; 
	my $al_minor_homo = ( $back{'hete2N_cnt_alleleTypeN'} > 1 ) ? $cnt_al{'homo_al'}{ $srt_al_homo[1] } : 0 ; 
	# $back{'hete2N_MinorAF'}           = ( $al_total_homo > 0 ) ? sprintf("%.4f", 100 * $al_minor_homo / $al_total_homo)  : -1 ; 
	$back{'hete2N_MinorAF'}           = ( $al_total_homo > 0 ) ? (100 * $al_minor_homo / $al_total_homo)  : -1 ; 
	$back{'hete2N_MinorAN'}           = $al_minor_homo; 
	$back{'hete2N_MajorAN'}           = ( $back{'hete2N_cnt_alleleTypeN'} > 0 ) ? $cnt_al{'homo_al'}{ $srt_al_homo[0] } : 0 ; 

	return(%back); 
}# cnt_type() 

# Example of rules : cnt_lmiss >= 295 
sub test_rule {
	my ($ah, $k, $c, $cut) = @_; 
	if ( $c eq '<' ) {
		return( $ah->{$k} < $cut ); 
	} elsif ( $c eq '<=' ) {
		return( $ah->{$k} <= $cut ); 
	} elsif ( $c eq '>' ) {
		return( $ah->{$k} > $cut ); 
	} elsif ( $c eq '>=') {
		return( $ah->{$k} >= $cut ); 
	} elsif ( $c eq '=' or $c eq '==' ) {
		return( $ah->{$k} == $cut ); 
	} else {
		die "Unknown [$c]\n"; 
	}
}# test_rule()




