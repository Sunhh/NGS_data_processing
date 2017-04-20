#!/usr/bin/perl
# 20170419 Generate separated alignments instead of concatenated one. 
use strict; 
use warnings; 
use LogInforSunhh; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"grp_list:s", 
	"pep_fas:s", "out_aln_pep:s", 
	"cds_fas:s", "out_aln_cds:s", 
	"tax_list:s", 
	"startColN:i", 
	"startRowN:i", 
	"min2max_var:f", 
	"minProt_len:i", 
	"para_trimal:s", "notTrimmed!", 
	"sepAln!", 
); 
$opts{'startColN'} //= 0; 
$opts{'startRowN'} //= 1; 
$opts{'min2max_var'} //= -1; 
$opts{'minProt_len'} //= -1; 
$opts{'para_trimal'} //= '-gt 0.9'; 


my $help_txt = <<HH; 

perl $0 -grp_list gene_grp_list -pep_fas all_pep.fas 

Required : muscle , trimal in the PATH

Format of -grp_list : gene1_tax1 \\t gene1_tax2 \\t gene1_tax3 .... 
Format of -pep_fas  : fasta file with name for 'gene1_tax1', 'gene1_tax2', ... 
Format of -cds_fas  : fasta file with name for 'gene1_tax1', 'gene1_tax2', ... 
Format of -tax_list : tax1 \\n tax2 \\n tax3 ... 

-tax_list       [''] If given, the conbined sequence has a better header than number. 
-startColN      [$opts{'startColN'}] 0-indexed. 
-startRowN      [$opts{'startRowN'}] 1-indexed. 
-min2max_var    [-1] [0-1]; 
-minProt_len    [-1] 

-para_trimal    ['-gt 0.9']
-notTrimmed     [Boolean] No trimming in trimal. 

-cds_fas        [filename]
-out_aln_cds    [filename] output .fasta file with cds. 

-out_aln_pep    [filename] Output .fasta file with protein sequences. This has been trimmed with trimal. 

-sepAln         [Boolean] Separate alignments by each gene. 

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
for (qw/grp_list pep_fas/) {
	defined $opts{$_} or &LogInforSunhh::usage($help_txt); 
}

my %in_seq = %{ $fs_obj->save_seq_to_hash( 'faFile'=>$opts{'pep_fas'}, 'has_head'=>1 ) }; 
&rmTailStar( \%in_seq ); 
my %in_cds; 
my $ofh_aln_cds; 
if ( defined $opts{'out_aln_cds'} ) {
	%in_cds = %{ $fs_obj->save_seq_to_hash( 'faFile'=>$opts{'cds_fas'}, 'has_head'=>1 ) }; 
	&rmTailStop( \%in_cds ); 
	$ofh_aln_cds = &openFH( $opts{'out_aln_cds'}, '>' ); 
}
my $ofh_aln_pep = \*STDOUT; 
if ( defined $opts{'out_aln_pep'}) {
	$ofh_aln_pep = &openFH( $opts{'out_aln_pep'}, '>' ); 
}

my @grp_list = @{ &load_gene_grp( $opts{'grp_list'}, $opts{'startColN'}, $opts{'startRowN'} ) }; 
my $cwd = &fileSunhh::_abs_path("."); 
my %taxName; 
defined $opts{'tax_list'} and %taxName = %{ &load_tax_list( $opts{'tax_list'} ) }; 

my %aln_seq; 
my %aln_cds; 
my %key_pep; 
my $used_cnt = 0; 
GRP:
for ( my $i=0; $i<@grp_list; $i++ ) {
	chdir($cwd); 
	my $tg = $grp_list[$i]; 
	my $tmp_dir = &fileSunhh::new_tmp_file(); 
	mkdir($tmp_dir); 
	&grpSeq_to_file( \%in_seq, $tg, "$tmp_dir/input.fa" ) and do { &tsmsg("[Wrn] Bad group skipped.\n"); &fileSunhh::_rmtree($tmp_dir); next GRP; }; 
	my ($min,$max) = &min_max_len_from_fas("$tmp_dir/input.fa"); 
	$min >= $opts{'minProt_len'} or do { &fileSunhh::_rmtree($tmp_dir); next GRP; }; 
	if ( $opts{'min2max_var'} >= 0 ) {
		$max <= $min * (1+$opts{'min2max_var'}) or do { &fileSunhh::_rmtree($tmp_dir); next GRP; }; 
	}
	chdir($tmp_dir); 
	&exeCmd_1cmd("muscle -in input.fa -out aln.fa") and do { &tsmsg("[Err] Skip bad group [$i]: @$tg\n"); chdir($cwd); &fileSunhh::_rmtree($tmp_dir); next GRP;}; 
	my %raw_aln_pep = %{ $fs_obj->save_seq_to_hash( 'faFile'=>'aln.fa', 'has_head'=>1 ) }; 
	for my $tk (keys %raw_aln_pep) { $raw_aln_pep{$tk}{'seq'} =~ s!\s!!g; } 


	my $str_cnum = ''; 
	$opts{'para_trimal'} =~ m!^\s*$! and $opts{'notTrimmed'} = 1; 
	if ( $opts{'notTrimmed'} ) {
		&exeCmd_1cmd("trimal -in aln.fa -out aln_trimal.fa $opts{'para_trimal'}"); 
	} else {
		$str_cnum = `trimal -in aln.fa -out aln_trimal.fa -colnumbering $opts{'para_trimal'}`; 
		chomp($str_cnum); $str_cnum =~ s!\s+!!g; 
	}
	my @kept_cnum = split( /,/, $str_cnum ); # These numbers are 0-based. 
	if ( $opts{'min2max_var'} >= 0 ) {
		my ($min2, $max2) = &min_max_len_from_fas("aln.fa"); 
		$max2 <= $min * (1+$opts{'min2max_var'}) or do { chdir($cwd); &fileSunhh::_rmtree($tmp_dir); next GRP; }; 
		my ($min3, $max3) = &min_max_len_from_fas("aln_trimal.fa"); 
		if ($min3 < $min) {
			$max2 <= $min3 * (1+$opts{'min2max_var'}) or do { chdir($cwd); &fileSunhh::_rmtree($tmp_dir); next GRP; }; 
		}
	}

	my $s2h = $fs_obj->save_seq_to_hash( 'faFile'=>"aln_trimal.fa", 'has_head'=>1 ); 
	chdir($cwd); 
	&fileSunhh::_rmtree($tmp_dir); 

	my %to_add_pep; 
	my %to_add_cds; 
	for (my $j=0; $j<@$tg; $j++) {
		my $seq_k = $tg->[$j]; 
		defined $s2h->{ $seq_k } or do { &tsmsg("[Err] Failed to find sequence for [$seq_k] in [$i] group: @$tg\n"); next GRP; }; 
		if ( defined $opts{'out_aln_cds'} ) {
			my $cdsByAA = &fastaSunhh::aa2cds_1seq( $in_cds{$seq_k}{'seq'}, $raw_aln_pep{$seq_k}{'seq'} ); 
			unless ( defined $cdsByAA ) {
				&tsmsg("[Err] Failed to match CDS sequence for [$seq_k] in [$i] group: @$tg\n"); 
				next GRP; 
			}
			$cdsByAA =~ s!\s!!g; 
			my $kept_cdsByAA = ''; 
			if ($opts{'notTrimmed'}) {
				$kept_cdsByAA = $cdsByAA; 
			} else {
				for my $ik (@kept_cnum) {
					$kept_cdsByAA .= substr($cdsByAA, $ik*3, 3); 
				}
			}
			$to_add_cds{$j} = $kept_cdsByAA; 
		}
		$s2h->{ $seq_k }{'seq'} =~ s/\s//gs; 
		$taxName{$j} //= "tax$j"; 
		$to_add_pep{$j} = $s2h->{ $seq_k }{'seq'}; 
	}
	$used_cnt ++; 
	for (my $j=0; $j<@$tg; $j++) {
		my $seq_k = $tg->[$j]; 
		push(@{$key_pep{$j}}, $seq_k); 
		$aln_seq{$j} .= "^^" . $to_add_pep{$j}; 
		if (defined $opts{'out_aln_cds'}) {
			$aln_cds{$j} .= "^^" . $to_add_cds{$j}; 
		}
	}
}

&tsmsg("[Rec] Total [$used_cnt] groups used.\n"); 

KN: 
for my $kn ( sort {$a<=>$b} keys %aln_seq ) {
	if ($opts{'sepAln'}) {
		$aln_seq{$kn} =~ s/\s//g; 
		$aln_seq{$kn} =~ s!^\^\^!!; 
		my @sub_pep = split(/\^\^/, $aln_seq{$kn}); 
		my @sub_cds ; 
		if ( defined $opts{'out_aln_cds'} ) {
			$aln_cds{$kn} =~ s!^\^\^!!; 
			@sub_cds = split(/\^\^/, $aln_cds{$kn});
		}
		for (my $i=0; $i<@sub_pep; $i++) {
			my $len = length( $sub_pep[$i] ); 
			$len > 0 or next KN; 
		}
		for (my $i=0; $i<@sub_pep; $i++) {
			my $len = length( $sub_pep[$i] ); 
			my $len_bar = ( $sub_pep[$i] =~ tr/-/-/ ); 
			my $len_aa  = $len - $len_bar; 
			$sub_pep[$i] =~ s!(.{60})!$1\n!g; chomp( $sub_pep[$i] ); 
			print {$ofh_aln_pep} ">$taxName{$kn}.$i $len aa [$len_aa]\n$sub_pep[$i]\n"; 
			if ( defined $opts{'out_aln_cds'} ) {
				$sub_cds[$i] =~ s!(.{60})!$1\n!g; chomp( $sub_cds[$i] ); 
				print {$ofh_aln_cds} ">$taxName{$kn}.$i\n$sub_cds[$i]\n"; 
			}
		}
		next KN; 
	}
	$aln_seq{$kn} =~ s/\s//g; 
	$aln_seq{$kn} =~ s!\^\^!!g; 
	my $len = length( $aln_seq{$kn} ); 
	my $len_bar = ( $aln_seq{$kn} =~ tr/-/-/ ); 
	my $len_aa  = $len-$len_bar; 
	$aln_seq{$kn} =~ s!(.{60})!$1\n!g; chomp( $aln_seq{$kn} ); 
	print {$ofh_aln_pep} ">$taxName{$kn} $len aa [$len_aa]\n$aln_seq{$kn}\n"; 
	if (defined $opts{'out_aln_cds'}) {
		$aln_cds{$kn} =~ s!\^\^!!g; 
		$aln_cds{$kn} =~ s!(.{60})!$1\n!g; chomp( $aln_cds{$kn} ); 
		print {$ofh_aln_cds} ">$taxName{$kn}\n$aln_cds{$kn}\n"; 
	}
}

sub load_gene_grp {
	# $_[0] : $opts{'grp_list'}
	# $_[1] : $opts{'startColN'}
	# $_[2] : $opts{'startRowN'}
	my @back; 
	my $fh = &openFH($_[0], '<'); 
	while (<$fh>) {
		$. >= $_[2] or next; 
		m/^\s*(#|$)/ and next; 
		chomp; 
		my @ta = split(/\t/, $_); 
		push(@back, [@ta[ $_[1] .. $#ta ]]); 
	}
	close($fh); 
	return( \@back ); 
}# load_gene_grp() 

sub grpSeq_to_file {
	# $_[0] : \%seq_hash 
	# $_[1] : \@id_list 
	# $_[2] : $out_file_name 
	my $is_bad = 0; 
	my $ofh = &openFH($_[2], '>'); 
	for my $tk (@{$_[1]}) {
		defined $_[0]->{$tk} or do { $is_bad = 1; last; }; 
		my $seq = $_[0]->{$tk}{'seq'}; 
		$seq =~ s/\s//gs; 
		$seq =~ s/(.{60})/$1\n/g; chomp($seq); 
		print {$ofh} ">$_[0]->{$tk}{'key'}\n$seq\n"; 
	}
	close($ofh); 
	if ($is_bad == 1) {
		return 1; 
	}
	return 0; 
}# grpSeq_to_file() 
sub load_tax_list {
	# $_[0] : tax_list
	my %back; 
	my $fh = &openFH($_[0], '<'); 
	my $cnt = 0; 
	while (<$fh>) {
		m/^\s*(#|$)/ and next; 
		chomp; 
		my @ta = split(/\t/, $_); 
		$back{$cnt} = $ta[0]; 
		$cnt ++; 
	}
	close($fh); 
	return \%back; 
}# load_tax_list() 
sub rmTailStar {
	# $_[0] : \%in_seq
	for (keys %{$_[0]}) {
		$_[0]{$_}{'seq'} =~ s/\s//g; 
		$_[0]{$_}{'seq'} =~ s/\*+$//; 
	}
	return; 
}
sub rmTailStop {
	for (keys %{$_[0]}) {
		$_[0]{$_}{'seq'} =~ s!\s!!g; 
		$_[0]{$_}{'seq'} = uc($_[0]{$_}{'seq'}); 
		my $ll = length( $_[0]{$_}{'seq'} ); 
		my $n_left = $ll % 3; 
		if ( $n_left == 1 ) {
			$_[0]{$_}{'seq'} .= 'NN'; 
		} elsif ( $n_left == 2 ) {
			$_[0]{$_}{'seq'} .= 'N'; 
		}
		$_[0]{$_}{'seq'} =~ s!(TAG|TAA|TGA)+$!!g; 
	}
	return; 
}
sub min_max_len_from_fas {
	# $_[0] : fasta file
	my ($min, $max); 
	my $s2h = $fs_obj->save_seq_to_hash( 'faFile'=>$_[0] ); 
	for my $tk (keys %$s2h) {
		$s2h->{$tk}{'seq'} =~ s/\s//gs; 
		my $l = length( $s2h->{$tk}{'seq'} );  
		$min //= $l; $min > $l and $min = $l; 
		$max //= $l; $max < $l and $max = $l; 
	}
	return($min, $max); 
}

