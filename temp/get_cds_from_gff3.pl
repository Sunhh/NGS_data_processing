#!/usr/bin/perl
### 20180504 Add [frame=\d] for futher usage. 
use strict; 
use warnings; 
use LogInforSunhh; 
use gffSunhh; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"genome_fas:s", 
	"genome_gff:s", 
	"help!", 
); 

my $help_txt = <<HH; 

perl $0 -genome_fas PG1All_v2_Scf.unmsk.fa -genome_gff r9_maker_good.gff3 > out.cds.fa

HH

( defined $opts{'genome_fas'} and defined $opts{'genome_gff'} ) or &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt);


my %str2num = qw(
 +     1
 -     -1
 1     1
 -1    -1
 0     -1
 plus  1
 minus -1
); 

my %cdsP_to_frame; 

my $gffF = $opts{'genome_gff'}; 
my $seqF = $opts{'genome_fas'}; 
my $oSeqWidth = 0; 
$oSeqWidth = int($oSeqWidth); 
my $log_detail = 0; 

&tsmsg("[Rec] Script begins.\n"); 

my $seqHR = $fs_obj->save_seq_to_hash( 'faFile'=>$seqF ); 
my %seq_hash; 
for my $tk (sort keys %$seqHR) {
	$seq_hash{$tk} = $seqHR->{$tk}{'seq'}; 
	$seq_hash{$tk} =~ s!\s!!g; 
}
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
		if ( keys %{$gff_hash{'lineN2hash'}{$offLnNum}{'attrib'}{'parentID'}} > 0 ) {
			my @ta1 = sort keys %{$gff_hash{'lineN2hash'}{$offLnNum}{'attrib'}{'parentID'}}; 
			my $ta1_txt = join(";", @ta1); 
			if ($top_name eq '') {
				$top_name = $ta1_txt; 
			} else {
				$top_name eq $ta1_txt or &stopErr("[Err] Unequal top_name : [ $top_name , $ta1_txt]\n"); 
			}
		}
		my @ta = split(/\t/, $gff_hash{'lineN2line'}{$offLnNum}); 
		my $cdsP_k = "$ta[3]:$ta[4]"; 
		$cdsP_to_frame{$cdsP_k} = $ta[7]; 
		$cdsP_to_frame{$cdsP_k} eq '.' and $cdsP_to_frame{$cdsP_k} = 0; 
		$cdsP_to_frame{$cdsP_k}++; 
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
	my $cdsP_1_k = "$posi_cds[0][0]:$posi_cds[0][1]"; 
	defined $cdsP_to_frame{$cdsP_1_k} or &stopErr("[Err] Failed to find frame for CDS position [$seq_hash{$top_chr}:$cdsP_1_k]\n"); 
	
	# get sequences. 
	my @sub_seqs; 
	for my $tr ( @posi_cds ) {
		push( @sub_seqs, substr($seq_hash{$top_chr}, $tr->[0]-1, $tr->[1]-$tr->[0]+1) ); 
	}
	$top_str == -1 and @sub_seqs = &rev_comp(@sub_seqs); 
	my $final_seq = join('', @sub_seqs); 
	$oSeqWidth > 0 and do { $final_seq =~ s!(.{$oSeqWidth})!$1\n!og; chomp($final_seq); }; 
	print STDOUT ">$top_name [frame=$cdsP_to_frame{$cdsP_1_k}]\n$final_seq\n"; 
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
