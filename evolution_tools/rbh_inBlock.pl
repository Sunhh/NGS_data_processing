#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 

	# Input: 
	"mcs_tbl:s", # ./ma_mo.collinearity.ks.tab
	"mcs_gff:s", # ma_mo.gff 
  "mcs_bp6:s@", # In file ma_mo.blast, there should not be ma_ma or mo_mo blast results. 

	"pl_rbh_byBp6:s", # /home/Sunhh/tools/github/NGS_data_processing/evolution_tools/rbh_byBp6.pl
  "para_rbh_byBp6:s", # ' -min_similarity 0.4 -max_lenDiffR 1.5 -log_lines 10000 '
); 

my $help_txt = <<HH; 

perl $0 -mcs_bp6 ma_mo.bp6 -mcs_bp6 mo_ma.bp6 -mcs_tbl ./ma_mo.collinearity.ks.tab -mcs_gff ./ma_mo.gff 

-pl_rbh_byBp6          [/home/Sunhh/tools/github/NGS_data_processing/evolution_tools/rbh_byBp6.pl]
-para_rbh_byBp6        [ -min_similarity 0.4 -max_lenDiffR 1.5 -log_lines 10000 ]

HH

&prepare_input(); 

sub prepare_input {
	$opts{'pl_rbh_byBp6'} //= '/home/Sunhh/tools/github/NGS_data_processing/evolution_tools/rbh_byBp6.pl'; 
	$opts{'para_rbh_byBp6'} //= ' -min_similarity 0.4 -max_lenDiffR 1.5 -log_lines 10000 '; 
	for my $tk (qw/mcs_tbl mcs_gff mcs_bp6/) {
		defined $opts{$tk} or do { &tsmsg("\n"); &tsmsg("[Err] Lack of parameter -$tk\n"); &LogInforSunhh::usage($help_txt); }; 
	}
	return; 
}# prepare_input() 

my @info_mcs_tbl = @{ &load_mcs_tbl($opts{'mcs_tbl'}) }; 
my %info_mcs_gff = %{ &load_mcs_gff($opts{'mcs_gff'}) }; 

my $wrk_dir = &fileSunhh::new_tmp_dir(); 
mkdir($wrk_dir) or &stopErr("[Err] Failed to create working dir [$wrk_dir]\n"); 
for my $t1 (@info_mcs_tbl) {
	my ( $blkID, $chr1,$chr1_S,$chr1_E, $chr2,$chr2_S,$chr2_E, $strand ) = @$t1; 
	&tsmsg("[Msg] Processing BlkID [$blkID][$chr1,$chr1_S,$chr1_E, $chr2,$chr2_S,$chr2_E]\n"); 

	(defined $info_mcs_gff{$chr1} and defined $info_mcs_gff{$chr2}) or &stopErr("[Err] $chr1 $chr2\n"); 
	my @gene1_list = map { $_->[0] } grep { !( $_->[1] > $chr1_E || $_->[2] < $chr1_S ) } @{$info_mcs_gff{$chr1}{'arr'}}; 
	my @gene2_list = map { $_->[0] } grep { !( $_->[1] > $chr2_E || $_->[2] < $chr2_S ) } @{$info_mcs_gff{$chr2}{'arr'}}; 
	&write_gene_list("$wrk_dir/list.g1", \@gene1_list); 
	&write_gene_list("$wrk_dir/list.g2", \@gene2_list); 
	my %gene1_hash = map { $_ => 1 } @gene1_list; 
	my %gene2_hash = map { $_ => 1 } @gene2_list; 

	&tsmsg("[Msg]   Generating $wrk_dir/new.bp6\n"); 
	my $ofh = &openFH("$wrk_dir/new.bp6", '>'); 
	for my $fn_bp6 (@{$opts{'mcs_bp6'}}) {
		&tsmsg("[Msg]   Processing bp6 [$fn_bp6]\n"); 
		my $fh = &openFH($fn_bp6, '<'); 
		my %cnt; 
		$cnt{'cntN_step'} = 100e3; 
		while (&wantLineC($fh)) {
			&fileSunhh::log_section($., \%cnt) and &tsmsg("[Msg]     $. line.\n"); 
			my @ta = &splitL("\t", $_); 
			( defined $gene1_hash{$ta[0]} and defined $gene2_hash{$ta[1]} ) or ( defined $gene1_hash{$ta[1]} and defined $gene2_hash{$ta[0]} ) or next; 
			print {$ofh} "$_\n"; 
		}
		close($fh); 
	}
	close($ofh); 

	# my $para_bp6 = join(' ', map { "-in_bp6 $_" } @{$opts{'mcs_bp6'}}); 
	my $para_bp6 = "-in_bp6 $wrk_dir/new.bp6"; 

	&exeCmd_1cmd("perl $opts{'pl_rbh_byBp6'} $para_bp6 -gene1_list $wrk_dir/list.g1 -gene2_list $wrk_dir/list.g2 -combine_gene12 $opts{'para_rbh_byBp6'} > $wrk_dir/blk_rbh.pairs"); 

	&tsmsg("[Msg]   Exporting $wrk_dir/blk_rbh.pairs\n"); 
	open F,'<',"$wrk_dir/blk_rbh.pairs" or &stopErr("[Err] Failed to open file [$wrk_dir/blk_rbh.pairs]\n"); 
	my @out_pairs; 
	while (&wantLineC(\*F)) {
		my @ta = &splitL("\t", $_); 
		(defined $info_mcs_gff{$chr1}{'hash'}{$ta[0]} and defined $info_mcs_gff{$chr2}{'hash'}{$ta[1]}) or @ta[0,1] = @ta[1,0]; 
		(defined $info_mcs_gff{$chr1}{'hash'}{$ta[0]} and defined $info_mcs_gff{$chr2}{'hash'}{$ta[1]}) or &stopErr("[Err] Bad pair: $ta[0] $ta[1]\n"); 
		push(@out_pairs, [ $blkID, $chr1, $ta[0], $chr2, $ta[1], $info_mcs_gff{$chr1}{'hash'}{$ta[0]}[0], $info_mcs_gff{$chr2}{'hash'}{$ta[1]}[0], $strand ]); 
	}
	close F; 
	for ( sort { $a->[5] <=> $b->[5] } @out_pairs ) {
		print STDOUT join("\t", @{$_}[0, 1,2,5, 3,4,6, 7])."\n"; 
	}
}

&fileSunhh::_rmtree($wrk_dir); 
&tsmsg("[Rec] $0 done.\n"); 


sub write_gene_list {
	my ($fn, $glist) = @_; 
	my $fh = &openFH($fn, '>'); 
	for (@$glist) {
		print {$fh} "$_\n"; 
	}
	close($fh); 
}# write_gene_list() 

# [Sunhh@wwz 02_ma_mo]$ deal_table.pl ./ma_mo.collinearity.ks.tab -col_head
# 0       BlkID
# 1       Chrom1
# 2       Start1
# 3       End1
# 4       Chrom2
# 5       Start2
# 6       End2
# 7       Strand
sub load_mcs_tbl {
	my ($fn) = @_; 
	my @back; 
	my $fh = &openFH($fn, '<'); 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		my ($blkID, $chr1,$chr1_S,$chr1_E, $chr2,$chr2_S,$chr2_E, $strand) = @ta[0, 1,2,3, 4,5,6, 7]; 
		$blkID eq 'BlkID' and next; 
		push(@back, [ $blkID, $chr1,$chr1_S,$chr1_E, $chr2,$chr2_S,$chr2_E, $strand ]); 
	}
	close($fh); 
	return(\@back); 
}# load_mcs_tbl() 

# [Sunhh@wwz 02_ma_mo]$ head -4 ma_mo.gff
# ma00    Cma_031679      55660   57727
# ma00    Cma_031673      57731   63801
sub load_mcs_gff {
	my ($fn) = @_; 
	my %back; 
	my $fh = &openFH($fn,'<'); 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		push(@{$back{$ta[0]}{'arr'}}, [@ta[1,2,3]]); 
		$back{$ta[0]}{'hash'}{$ta[1]} = [@ta[2,3]]; 
	}
	close($fh); 
	for (keys %back) {
		@{$back{$_}{'arr'}} = sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @{$back{$_}{'arr'}}; 
	}
	return(\%back); 
}# load_mcs_gff() 


