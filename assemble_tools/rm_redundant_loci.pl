#!/usr/bin/perl
# 20150121 Edit bug when the alignment block start with 1
use strict; 
use warnings; 
use LogInforSunhh; 
use fastaSunhh; 
my $fs_obj = fastaSunhh->new(); 
use fileSunhh; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"evalue:s", "minLen:i", "num_threads:i", "similarity:f", 
	"dbFas:s", "inFas:s", "outFas:s", 
	"help!"
); 

$opts{'evalue'} //= 1e-5; 
#$opts{'evalue'} = $opts{'evalue'} + 0; 
$opts{'num_threads'} //= 1; 
$opts{'minLen'} //= 50; 
$opts{'similarity'} //= 0; 

$opts{'blastn'} //= 'blastn'; 
$opts{'mkdb'}   //= 'makeblastdb'; 
$opts{'blastn_para'} //= '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand" -task blastn'; 

my $help_txt = <<HH; 

perl $0   -inFas input_nucl.fasta   -outFas output_nucl_rmR.fasta

-help

-evalue        [$opts{'evalue'}]
-num_threads   [$opts{'num_threads'}]
-minLen        [$opts{'minLen'}]
-similarity    [$opts{'similarity'}] 

-blastn_para   ['$opts{'blastn_para'}']

HH


defined $opts{'inFas'} or &LogInforSunhh::usage($help_txt); 
defined $opts{'outFas'} or &LogInforSunhh::usage($help_txt); 
defined $opts{'help'} and &LogInforSunhh::usage($help_txt); 

my %qryFa_h = %{ $fs_obj->save_seq_to_hash( 'faFile' => $opts{'inFas'} ) }; 
my %dbFa_h = (); 
defined $opts{'dbFas'} and %dbFa_h = %{ $fs_obj->save_seq_to_hash( 'faFile' => $opts{'dbFas'} ) }; 
my $oFh = &openFH($opts{'outFas'}, '>'); 


while ( keys %qryFa_h > 0 ) {
	&rm_redundant(\%qryFa_h, \%dbFa_h); 
	my $n = keys %dbFa_h; 
}

for my $tk (sort { $dbFa_h{$a}{'Order'} <=> $dbFa_h{$b}{'Order'} } keys %dbFa_h) {
	chomp( $dbFa_h{$tk}{'seq'} ); 
	print {$oFh} ">$dbFa_h{$tk}{'head'}\n$dbFa_h{$tk}{'seq'}\n"; 
}

sub rm_redundant {
	# $_[0] : \%qryFa_h 
	# $_[1] : \%dbFa_h
	$_[1] //= {}; 

	if (keys %{$_[0]} == 0) {
		&tsmsg("[Wrn] No qryFa needs to be treated.\n"); 
		return ($_[1]); 
	}
	my $tmpF_sbj = fileSunhh::new_tmp_file(); 
	if (keys %{$_[1]} == 0) {
		my ($lk) = &write_longest_to_file( $_[0], $tmpF_sbj ); 
		$_[1]->{$lk} = $_[0]->{$lk}; 
		delete $_[0]->{$lk}; 
	} else {
		&write_longest_to_file( {}, $tmpF_sbj, $_[1] ); 
	}

	my $tmpF_qry = fileSunhh::new_tmp_file(); 
	if ( keys %{$_[0]} == 0 ) {
		# Finished. 
		return ($_[1]); 
	}
	my ($qk) = &write_longest_to_file( $_[0], $tmpF_qry ); 
	&tsmsg("[Msg] Treating [$qk]\n"); 
	# delete $_[0]->{$qk} ; 
	my $tmpF_bn6 = "${tmpF_sbj}.bn6"; 
	&run_blastn( $tmpF_sbj, $tmpF_qry, $tmpF_bn6 ); 
	my %q_Rloc = %{ &r_loc_in_bn6( $tmpF_bn6 ) }; 
	$q_Rloc{$qk} //= []; 
	if ( @{$q_Rloc{$qk}} == 0 ) {
		my $ss = $_[0]->{$qk}{'seq'}; chomp($ss); 
		my $text = ">$_[0]->{$qk}{'head'}\n$ss\n"; 
		&fileSunhh::write2file( $tmpF_sbj, $text, '>>' ); 
	} else {
		my $qlen = $_[0]->{$qk}{'len'}; 
		my @q_Uloc = @{ &rest_loc( $q_Rloc{$qk}, $qlen ) }; 
		my $q_seq = $_[0]->{$qk}{'seq'}; $q_seq =~ s/\s//g; 
		for my $tr1 (@q_Uloc) {
			my ($ts, $te) = ($tr1->[0], $tr1->[1]); 
			$te-$ts+1 >= $opts{'minLen'} or next; 
			my $ss = substr( $q_seq, $ts-1, $te-$ts+1 ); 
			$ss =~ s/(.{70})/$1\n/g; chomp($ss); 
			my $text = ">${qk}__${ts}_${te}\n$ss\n"; 
			&fileSunhh::write2file( $tmpF_sbj, $text, '>>' ); 
		}
	}

	%{$_[1]} = %{$fs_obj->save_seq_to_hash( 'faFile' => $tmpF_sbj ) }; 
	delete($_[0]->{$qk}); 

	&exeCmd_1cmd("rm -f $tmpF_sbj $tmpF_qry $tmpF_bn6"); 
	return ($_[1]); 
}# rm_redundant

sub rest_loc {
	# $_[0] : [ [s1, e1], [s2, e2], ... ]
	# $_[1] : length >= max(e)
	unless (defined $_[1] and $_[1] > 0) {
		for my $tr (@{$_[0]}) {
			$_[1] //= $tr->[1]; 
			$_[1] < $tr->[1] and $_[1] = $tr->[1]; 
		}
	}
	$_[1] > 0 or &stopErr("[Err] No good length.\n"); 

	@{$_[0]} == 0 and return [[1, $_[1]]]; 
	my $prevE = 0; 
	@{$_[0]} = sort { $a->[0] <=> $b->[0] } @{$_[0]}; 
	my @back; 
	for my $tr1 (@{$_[0]}) {
		if ($tr1->[0] > $prevE+1) {
			push(@back, [$prevE+1, $tr1->[0]-1]); 
			$prevE = $tr1->[1]; 
		}
		$prevE < $tr1->[1] and $prevE = $tr1->[1]; 
	}
	$prevE < $_[1] and push(@back, [$prevE+1, $_[1]]); 
	return(\@back); 
}

sub r_loc_in_bn6 {
	# $_[0] : in.bn6 
	my $fh = &openFH($_[0], '<'); 
	my %sep_blk; 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		$ta[2] >= $opts{'similarity'} or next; 
		$ta[9] < $ta[8] and @ta[8,9]=@ta[9,8]; 
		$ta[9]-$ta[8]+1 >= $opts{'minLen'} or next; 
		push(@{$sep_blk{$ta[1]}}, [@ta[8,9]]); 
	}
	close($fh); 
	for my $tr (keys %sep_blk) {
		my $tb = $ms_obj->mergeLocBlk( $sep_blk{$tr} , 'dist2join'=>1 ); 
		@{$sep_blk{$tr}} = @$tb; 
	}
	return (\%sep_blk); 
}# r_loc_in_bn6() 

sub q_loc_in_bn6 {
	# $_[0] : in.bn6 
	my $fh = &openFH($_[0], '<'); 
	my %sep_blk; 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		$ta[2] >= $opts{'similarity'} or next; 
		$ta[6] > $ta[7] and @ta[6,7]=@ta[7,6]; 
		$ta[7]-$ta[6]+1 >= $opts{'minLen'} or next; 
		push(@{$sep_blk{$ta[0]}}, [@ta[6,7]]); 
	}
	close($fh); 
	for my $tr (keys %sep_blk) {
		my $tb = $ms_obj->mergeLocBlk( $tr , 'dist2join'=>1 ); 
		@$tr = @$tb; 
	}
	return (\%sep_blk); 
}# q_loc_in_bn6() 

sub run_blastn {
	# $_[0] : qry_file 
	# $_[1] : sbj_file
	# $_[2] : out_bn6 
	$_[2] //= "$_[0].bn6"; 
	&exeCmd_1cmd( "$opts{'mkdb'} -dbtype nucl -in $_[1]" ); 
	&exeCmd_1cmd( "$opts{'blastn'} $opts{'blastn_para'} -evalue $opts{'evalue'} -num_threads $opts{'num_threads'} -query $_[0] -db $_[1] -out $_[2]" ); 
	&exeCmd_1cmd( "rm -f $_[1]*.nhr $_[1]*.nin $_[1]*.nsq $_[1]*.nal" ); 
	return (); 
}# run_blastn 

sub write_longest_to_file {
	# $_[0] : \%qryFa_h 
	# $_[1] : $out_file_name 
	# $_[2] : \%addFa_h \
	$_[2] //= {}; 
	if ( (keys %{$_[0]}) == 0 and (keys %{$_[2]} == 0) ) {
		&stopErr("[Err] Both input are empty.\n"); 
	}
	my $ofh = &openFH($_[1], '>'); 
	if ( (keys %{$_[2]}) > 0 ) {
		for my $ak (sort { $_[2]->{$a} <=> $_[2]->{$b} } keys %{$_[2]}) {
			print {$ofh} ">$_[2]->{$ak}{'head'}\n$_[2]->{$ak}{'seq'}\n"; 
		}
	}

	if ( (keys %{$_[0]}) == 0 ) {
		close($ofh); 
		return undef(); 
	}
	for my $tk (keys %{$_[0]}) {
		unless ( defined $_[0]{$tk}{'len'} ) {
			my $ss = $_[0]{$tk}{'seq'}; 
			$ss =~ s/\s//g; 
			$_[0]{$tk}{'len'} = length($ss); 
		}
	}
	my ($lk) = sort { $_[0]->{$b}{'len'} <=> $_[0]->{$a}{'len'} || $_[0]->{$a}{'Order'} <=> $_[0]->{$a}{'Order'} } keys %{$_[0]}; 

	print {$ofh} ">$_[0]->{$lk}{'head'}\n$_[0]->{$lk}{'seq'}\n"; 
	close($ofh); 
	return ($lk); 
}# write_longest_to_file ()

