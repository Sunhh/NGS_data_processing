#!/usr/bin/perl
# CdEx: Coding exon part. 
# I want to remove these short coding exon part from gff3 files to get clean and good prot2genome alignments. 
# Remove from both ends. 
use strict; 
use warnings; 
use LogInforSunhh; 
use mathSunhh; 
my $ms = mathSunhh->new(); 

-t and !@ARGV and die "perl $0 CM_CS_WM97_Arab_SprotPln.spaln_M4.fmt.gff3\n"; 

my $shrtCdEx_len = 30; 
my $singleEx_dist = 0; 

my @geneLines; 

while (<>) {
	chomp; 
	if (m!^\s*(#|$)!) {
		print STDOUT "$_\n"; 
		next; 
	}
	my @ta = split(/\t/, $_); 
	if ( $ta[2] =~ m/^protein_match$/i ) {
		if (@geneLines > 0) {
			&rmShrt(\@geneLines, $shrtCdEx_len, $singleEx_dist, 'left'); 
			&rmShrt(\@geneLines, $shrtCdEx_len, $singleEx_dist, 'right'); 
			&outGff(\@geneLines, \*STDOUT); 
		}
		@geneLines = (); 
		push(@geneLines, [@ta]); 
	} elsif ( $ta[2] =~ m/^match_part$/i ) {
		push(@geneLines, [@ta]); 
	} else {
		die "$_\n"; 
	}
}
if (@geneLines > 0) {
	&rmShrt(\@geneLines, $shrtCdEx_len, $singleEx_dist, 'left'); 
	&rmShrt(\@geneLines, $shrtCdEx_len, $singleEx_dist, 'right'); 
	&outGff(\@geneLines, \*STDOUT); 
}
@geneLines = (); 

sub outGff ($$) {
	my ($gl_aref, $ofh) = @_; 
	$ofh //= \*STDOUT; 
	for (my $i=0; $i<@$gl_aref; $i++) {
		print {$ofh} join("\t", @{$gl_aref->[$i]})."\n"; 
	}
	return; 
}

sub rmShrt($$$$$$) {
	my ($gl_aref, $len, $dist, $left_right, $chk_type_aref, $top_type_aref) = @_; 
	$len //= 30; 
	$dist //= 2000; 
	$left_right //= 'left'; 
	$chk_type_aref //= [qw/match_part/]; 
	$top_type_aref //= [qw/protein_match/]; 
	my %is_chk; 
	my %is_top; 
	$is_chk{$_} = 1 for (@$chk_type_aref) ; 
	$is_top{$_} = 1 for (@$top_type_aref) ; 
	
	# Check if need to be trimmmed. 
	my %rm_idx; 
	my @ex_len; # (Ex_len)
	my @ex_info; # [ExS, ExE, geneLine_idx]
	for ( my $i=0; $i<@$gl_aref; $i++ ) {
		defined $is_chk{$gl_aref->[$i][2]} or next; 
		push(@ex_len, $gl_aref->[$i][4]-$gl_aref->[$i][3]+1); 
		push( @ex_info, [$gl_aref->[$i][3], $gl_aref->[$i][4], $i] ); 
	}
	
	if ( $left_right eq 'left' ) {
	} elsif ( $left_right eq 'right' ) {
		@ex_len = reverse(@ex_len); 
		@ex_info = reverse(@ex_info); 
	} else {
		&stopErr("[Err]unknown left_right=[$left_right]\n"); 
	}
	if ( @ex_len == 0 or ( @ex_len == 1 and $ex_len[0] <= $len )) {
		@$gl_aref = (); 
		return; 
	}
	
	if ($ex_len[0] <= $len) {
		defined $ex_info[1][0] or die "@{$ex_info[0]}\n"; 
		my $intron_len = &intron_len( $ex_info[0][0], $ex_info[0][1], $ex_info[1][0], $ex_info[1][1] ); 
		$intron_len >= $dist and $rm_idx{$ex_info[0][2]} = 1; 
	}
	for (my $i=1; $i<@ex_len; $i++) {
		my $stat = $ms->ins_calc( [ @ex_len[0..$i] ] ); 
		$stat->{'MEAN'} > $len and last; 
		$rm_idx{$ex_info[$i][2]} = 1; 
	}
	# Do not remove exon longer than $len in the other end. 
	for (my $i=$#ex_len; $i>=0; $i--) {
		defined $rm_idx{$ex_info[$i][2]} or next; 
		$rm_idx{$ex_info[$i][2]} == 1 or next; 
		$ex_len[$i] <= $len and last; 
		$ex_len[$i] > $len and $rm_idx{$ex_info[$i][2]} = 0; 
	}
	
	# Now remove short exons. 
	my @filtered_gl; 
	my ($gl_min, $gl_max); 
	for (my $i=0; $i<@$gl_aref; $i++) {
		defined $rm_idx{$i} and $rm_idx{$i} == 1 and next; 
		push(@filtered_gl, $gl_aref->[$i]); 
		$gl_min //= $gl_aref->[$i][3]; 
		$gl_min = $ms->min($gl_min, $gl_aref->[$i][3]); 
		$gl_max //= $gl_aref->[$i][4]; 
		$gl_max = $ms->max($gl_max, $gl_aref->[$i][4]); 
	}
	# Give right boundary of 'protein_match'
	for (my $i=0; $i<@filtered_gl; $i++) {
		defined $is_top{$filtered_gl[$i][2]} or next; 
		$filtered_gl[$i][3] < $gl_min and $filtered_gl[$i][3] = $gl_min; 
		$filtered_gl[$i][4] > $gl_max and $filtered_gl[$i][4] = $gl_max; 
	}
	
	@$gl_aref = @filtered_gl; 
	return; 
}# rmShrt() 

sub intron_len {
	my ($s1,$e1,$s2,$e2) = @_; 
	return ($ms->min(abs($s1-$s2), abs($s1-$e2), abs($e1-$s2), abs($e1-$e2)) + 1); 
}

