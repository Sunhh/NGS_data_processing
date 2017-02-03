#!/usr/bin/perl
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
	"pep_fas:s", 
	"tax_list:s", 
	"startColN:i", 
	"startRowN:i", 
	"min2max_var:f", 
	"minProt_len:i", 
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
Format of -tax_list : tax1 \\n tax2 \\n tax3 ... 

-tax_list       [''] If given, the conbined sequence has a better header than number. 
-startColN      [$opts{'startColN'}] 0-indexed. 
-startRowN      [$opts{'startRowN'}] 1-indexed. 
-min2max_var    [-1] [0-1]; 
-minProt_len    [-1] 

-para_trimal    ['-gt 0.9']

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
for (qw/grp_list pep_fas/) {
	defined $opts{$_} or &LogInforSunhh::usage($help_txt); 
}

my %in_seq = %{ $fs_obj->save_seq_to_hash( 'faFile'=>$opts{'pep_fas'}, 'has_head'=>1 ) }; 
&rmTailStar( \%in_seq ); 
my @grp_list = @{ &load_gene_grp( $opts{'grp_list'}, $opts{'startColN'}, $opts{'startRowN'} ) }; 
my $cwd = &fileSunhh::_abs_path("."); 
my %taxName; 
defined $opts{'tax_list'} and %taxName = %{ &load_tax_list( $opts{'tax_list'} ) }; 

my %aln_seq; 
my $used_cnt = 0; 
for ( my $i=0; $i<@grp_list; $i++ ) {
	chdir($cwd); 
	my $tg = $grp_list[$i]; 
	my $tmp_dir = &fileSunhh::new_tmp_file(); 
	mkdir($tmp_dir); 
	&grpSeq_to_file( \%in_seq, $tg, "$tmp_dir/input.fa" ) and do { &tsmsg("[Wrn] Bad group skipped.\n"); &fileSunhh::_rmtree($tmp_dir); next; }; 
	my ($min,$max) = &min_max_len_from_fas("$tmp_dir/input.fa"); 
	$min >= $opts{'minProt_len'} or do { &fileSunhh::_rmtree($tmp_dir); next; }; 
	if ( $opts{'min2max_var'} >= 0 ) {
		$max <= $min * (1+$opts{'min2max_var'}) or do { &fileSunhh::_rmtree($tmp_dir); next; }; 
	}
	chdir($tmp_dir); 
	&exeCmd_1cmd("muscle -in input.fa -out aln.fa"); 
	&exeCmd_1cmd("trimal -in aln.fa -out aln_trimal.fa $opts{'para_trimal'}"); 
	if ( $opts{'min2max_var'} >= 0 ) {
		my ($min2, $max2) = &min_max_len_from_fas("aln.fa"); 
		$max2 <= $min * (1+$opts{'min2max_var'}) or do { chdir($cwd); &fileSunhh::_rmtree($tmp_dir); next; }; 
		my ($min3, $max3) = &min_max_len_from_fas("aln_trimal.fa"); 
		if ($min3 < $min) {
			$max2 <= $min3 * (1+$opts{'min2max_var'}) or do { chdir($cwd); &fileSunhh::_rmtree($tmp_dir); next; }; 
		}
	}

	my $s2h = $fs_obj->save_seq_to_hash( 'faFile'=>"aln_trimal.fa", 'has_head'=>1 ); 
	chdir($cwd); 
	&fileSunhh::_rmtree($tmp_dir); 

	$used_cnt ++; 
	for (my $j=0; $j<@$tg; $j++) {
		defined $s2h->{ $tg->[$j] } or &stopErr("[Err] Failed to find sequence for [$tg->[$j]]\n"); 
		$s2h->{ $tg->[$j] }{'seq'} =~ s/\s//gs; 
		$taxName{$j} //= "tax$j"; 
		$aln_seq{$j} .= $s2h->{ $tg->[$j] }{'seq'}; 
	}
}

&tsmsg("[Rec] Total [$used_cnt] groups used.\n"); 

for my $kn ( sort {$a<=>$b} keys %aln_seq ) {
	$aln_seq{$kn} =~ s/\s//g; 
	my $len = length( $aln_seq{$kn} ); 
	my $len_bar = ( $aln_seq{$kn} =~ tr/-/-/ ); 
	my $len_aa  = $len-$len_bar; 
	$aln_seq{$kn} =~ s!(.{60})!$1\n!g; chomp( $aln_seq{$kn} ); 
	print STDOUT ">$taxName{$kn} $len aa [$len_aa]\n$aln_seq{$kn}\n"; 
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

