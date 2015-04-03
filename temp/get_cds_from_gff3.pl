#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use gffSunhh; 

!@ARGV and die "perl $0 fasdfa\n"; 

my %str2num = qw(
 +     1
 -     -1
 1     1
 -1    -1
 0     -1
 plus  1
 minus -1
); 

my $gffF = 'r9_maker_good.gff3'; 
my $seqF = 'PG1All_v2_Scf.unmsk.fa'; 
my $oSeqWidth = 0; 
$oSeqWidth = int($oSeqWidth); 
my $log_detail = 0; 

&tsmsg("[Rec] Script begins.\n"); 

open FA,'<',"$seqF" or die; 
my %seq_hash; 
{
my ($tk)=(''); 
while (<FA>) {
	chomp; 
	if (m/^\s*\>(\S+)/) {
		$tk = $1; 
	}else{
		$seq_hash{$tk} .= $_; 
	}
}
$seq_hash{$_} =~ s!!!g foreach (keys %seq_hash); 
}
close FA; 
&tsmsg("[Rec] seqfile [$seqF] read in.\n"); 

my $gs = gffSunhh->new(); 
my ($gff_href, $seq_href) = $gs->read_gff3File('gffFile'=>$gffF, 'saveFa'=>0, 'top_hier'=>{ 'mrna'=>1, 'match'=>2, 'protein_match'=>3, 'expressed_sequence_match'=>4 }); 
&tsmsg("[Rec] gfffile [$gffF] read in.\n"); 

my %gff_hash = %$gff_href; undef( $gff_href ); 
TOPID: 
for my $topID ( keys %{$gff_hash{'lineN_group'}} ) {
	# &tsmsg("[Msg] Dealing topID [$topID]\n"); 
	my @posi_top; 
	my @posi_cds; 
	my $curLnNum = $gff_hash{'lineN_group'}{$topID}{'curLn'}[0]; 
	my $top_str = ''; 
	my $top_chr = ''; 
	my $top_name = ''; 
	if ( defined $gff_hash{'lineN2hash'}{$curLnNum} ) {
		$gff_hash{'lineN2hash'}{$curLnNum}{'type'} =~ m/^mRNA$/i or next TOPID; 
		$top_str = lc( $gff_hash{'lineN2hash'}{$curLnNum}{'strand'} ); 
		defined $str2num{$top_str} or $top_str = ''; 
		$top_str = $str2num{$top_str} // ''; 
		$top_chr = $gff_hash{'lineN2hash'}{$curLnNum}{'seqID'} // ''; 
		$top_name = $gff_hash{'lineN2hash'}{$curLnNum}{'attrib'}{'featID'} // ''; 
		@posi_top = ([ $gff_hash{'lineN2hash'}{$curLnNum}{'start'}, $gff_hash{'lineN2hash'}{$curLnNum}{'end'} ]); 
	}
	OFFID: 
	for my $offLnNum ( @{$gff_hash{'lineN_group'}{$topID}{'offLn'}} ) {
		defined $gff_hash{'lineN2hash'}{$offLnNum} or next OFFID; 
		$gff_hash{'lineN2hash'}{$offLnNum}{'type'} =~ m/^CDS$/i or next OFFID; 
		$top_str eq '' and $top_str = $str2num{ lc($gff_hash{'lineN2hash'}{$offLnNum}{'strand'}) } // ''; 
		$top_chr eq '' and $top_chr = $gff_hash{'lineN2hash'}{$offLnNum}{'seqID'} // ''; 
		$top_name eq '' and $top_name = $gff_hash{'lineN2hash'}{$offLnNum}{'attrib'}{'parentID'} // ''; 
		push( @posi_cds, [ $gff_hash{'lineN2hash'}{$offLnNum}{'start'}, $gff_hash{'lineN2hash'}{$offLnNum}{'end'} ] ); 
	}
	$top_str eq '' and do { &tsmsg("[Wrn] No strand information for topID=[$topID]\n"); $top_str = 1; }; 
	$top_chr eq '' and do { &stopErr("[Err] No top_chr found for [$topID]\n"); };
	defined $seq_hash{$top_chr} or do { &tsmsg( "[Wrn] No [$top_chr] sequence found for [$topID]\n" ); next TOPID; }; 
	# Check @posi_cds 
	for my $tr (@posi_cds) {
		$tr->[0] > $tr->[1] and &stopErr("[Err] I can't accept start[$tr->[0]] > end[$tr->[1]] in gff3 file.\n"); 
	}
	if ($top_str == -1) {
		@posi_cds = sort { $b->[0] <=> $a->[0] } @posi_cds; 
	}
	
	# get sequences. 
	my @sub_seqs; 
	for my $tr ( @posi_cds ) {
		push( @sub_seqs, substr($seq_hash{$top_chr}, $tr->[0]-1, $tr->[1]-$tr->[0]+1) ); 
	}
	$top_str == -1 and @sub_seqs = &rev_comp(@sub_seqs); 
	my $final_seq = join('', @sub_seqs); 
	$oSeqWidth > 0 and do { $final_seq =~ s!(.{$oSeqWidth})!$1\n!og; chomp($final_seq); }; 
	print STDOUT ">$top_name\n$final_seq\n"; 
}

sub rev_comp {
	my @back; 
	for (@_) {
		my $ts = reverse($_); 
		$ts =~ tr/acgturykmbvdhACGTURYKMBVDHwWsSnN/tgcaayrmkvbhdTGCAAYRMKVBHDwWsSnN/; 
		push(@back, $ts); 
	}
	return @back; 
}
