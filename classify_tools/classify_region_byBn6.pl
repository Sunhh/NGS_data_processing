#!/usr/bin/env perl
# By Sunhh. (hs738@cornell.edu)
# Use blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand staxids sscinames sskingdoms stitle'
#   result to type regions in the query into Eukaryota or not.
# 2013-10-17 Version 1.
# 2013-11-01 Edit to assign Include/Exclude kingdoms.
# 2014-03-04 Fix a bug in which we fail to classify some end-to-end alignments. 
# [4/21/2022] Join blocks according to InEx classes instead of kingdom classes.
use strict;
use warnings;
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long;
my %opts;

GetOptions(
	\%opts,
	"byHsp!",
	"joinInEx:s", "maxUn:i",
	"out:s",
	"InPlastid!",
	"InList:s", "ExList:s",
	"txid2KingList:s",
	"help!"
);

sub writeFH {
	my $fname = shift;
	my $tfh;
	open $tfh,'>',"$fname" or die "$!\n";
	return $tfh;
}
# oFmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand staxids sscinames sskingdoms stitle'
#         0      1      2      3      4        5       6      7    8      9    10     11       12   13   14      15      16        17         18
# Pre configuration.
my $only1Hit = ($opts{byHsp}) ? 0 : 1 ; # 如果指定这个参数为1 (Not using -byHsp parameter), 同一个query区间, 只会计算一次同一个Hit的cover, 这样如果query在Hit内重复出现(多个hsp), 这个Hit对该query单元的支持贡献也只有一次; 关闭为1; Not used now. 
$opts{maxUn} //= 1;

my %skingdom=qw(
NA        0
Viruses   1
Bacteria  2
Archaea   3
rDNA           4
Chloroplast    5
Mitochondrion  6
Plastid        7
Eukaryota      8
Satellite      9
Bacteria;Eukaryota  10
Eukaryota;Viruses   11
In                  10000
Ex                  10001
);
my %isInEx=qw(
NA        Ex
Viruses   Ex
Bacteria  Ex
Archaea   Ex
Eukaryota In
Satellite In
Chloroplast    Ex
Mitochondrion  Ex
Plastid        Ex
rDNA           Ex
); # Should be excluded or included.

my %txid2King; 
if (defined $opts{txid2KingList}) {
	my $tfh = &openFH($opts{txid2KingList}, '<'); 
	while (<$tfh>) {
		my ($id,$kn) = (split)[0,1]; 
		defined $skingdom{$kn} or $skingdom{$kn} = -100; 
		$txid2King{$id} = $kn; 
	}
	close $tfh; 
}#End if -txid2KingList


if ($opts{InPlastid}) {
	$isInEx{'Chloroplast'} = 'In';
	$isInEx{'Mitochondrion'} = 'In';
	$isInEx{'Plastid'} = 'In';
}
if (defined $opts{InList}) {
	for my $tK (split(/:/, $opts{InList})) {
		defined $skingdom{$tK} or do { &tsmsg("[Wrn] No kingdom [$tK] defined.\n"); };
		$isInEx{$tK} = 'In';
		&tsmsg("[Rec]Subject kingdom [$tK] is included.\n");
	}
}
if (defined $opts{ExList}) {
	for my $tK (split(/:/, $opts{ExList})) {
		defined $skingdom{$tK} or do { &tsmsg("[Wrn]No kingdom [$tK] defined.\n"); };
		$isInEx{$tK} = 'Ex';
		&tsmsg("[Rec]Subject kingdom [$tK] is excluded.\n");
	}
}


sub usage {
	my $version = "v1.0";
	my $last_time = "2013-10-17";
	$last_time = '2013-10-24';
	my $tKK = join( ":", sort { $skingdom{$a} <=> $skingdom{$b} } keys %skingdom );
	my $tII = join( ":", sort { $skingdom{$a} <=> $skingdom{$b} } grep { $isInEx{$_} eq "In" } keys %isInEx );
	my $tEE = join( ":", sort { $skingdom{$a} <=> $skingdom{$b} } grep { $isInEx{$_} eq "Ex" } keys %isInEx );

	print STDOUT <<INFO;
##########################################################################
# perl $0 in.bn6
##########################################################################
# Type query regions by the blastn to nt result.
# oFmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand staxids sscinames sskingdoms stitle'
##########################################################################
# Self use currently.
# Version       : $version
# Last modified : $last_time
# Contact       : hs738\@cornell.edu
##########################################################################
# -help
#
# -byHsp     Query region coverage is calculated by HSP instead of by Hit.
#             The default mode is by Hit.
#             Not used now.
#
# -txid2KingList [filename] File to tell sbjct txid to kingdom class. Format: TXID\\tKingdom\n
#
# -joinInEx  [out file name] joined adjacent In/Excluding units with similar types.
# -maxUn     [integer][$opts{maxUn}] Max number of "Un" type gap between adjacent typed regions.
# -InPlastid [Boolean] Include plastids instead of excluding them.
#             Could be replaced by -InList & -ExList;
#
# -InList    [String] [Kingdom1:Kingdom2:...] Kingdoms to be included.
# -ExList    [String] [Kingdom1:Kingdom2:...] Kingdoms to be excluded. Overwrite -InList
##########################################################################
# Kingdom list: $tKK
# Included    : $tII
# Excluded    : $tEE
#   NA here means the query hits a database sequence with no taxonomy information.
#   Usually they are artificial sequences.
##########################################################################
INFO
	exit 1;
}#End for usage.

if (-t and !@ARGV) {
	&usage();
} elsif ( $opts{help} ) {
	&usage();
}

my $is_joinInEx=0;
my $joinInExFH;
defined $opts{joinInEx} and do { $joinInExFH = &writeFH( $opts{joinInEx} ); $is_joinInEx = 1; };
my $oFH = \*STDOUT;
defined $opts{out} and $oFH = &writeFH( $opts{out} );


# oFmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand staxids sscinames sskingdoms stitle'
#         0      1      2      3      4        5       6      7    8      9    10     11       12   13   14      15      16        17         18
# Pre configuration.

my %qBlks; # {qid}=>[qstart, qend, skingdom_class]
my %qLen; # {qid}=>query_length
my %dvd_site; # {qid}{block-edge}=>1
my %scf_order; 
# Divide query into minimum block units that should be intact in any alignment.
while (<>) {
	/^qseqid\t/ and next;
	chomp;
	my @ta = split(/\t/, $_);
	#### Do some filter if needed to remove un-reliable alignments.
	# Filtering work should be done here instead of later.
	# I prefer to use all acceptable hsp instead of top hits or top hsp, because the sequence length varies and because I limit number of hits when blasting.
	#### Record query blocks assigned to a kingdom_class.
	my $qid = $ta[0]; # Query name
	$scf_order{$qid} //= $.; 
	my $sid = $ta[1]; # Hit name
	my $skd = $ta[17]; # Hit kingdom class
	my $sTXID = $ta[15]; 
	defined $txid2King{$sTXID} and $skd = $txid2King{$sTXID}; 
	my $stitle = $ta[18]; # Hit definition line.
	if      ( $skd eq 'N/A' ) {
		$skd='NA';
	} elsif ( $stitle =~ m/\bchloroplast\b/i ) {
		# $stitle !~ m/ribosom/i and $skd='Chloroplast';
		$skd='Chloroplast';
	} elsif ( $stitle =~ m/\b(?:mitochondrial|mitochondrion)\b/i ) {
		# $stitle !~ m/ribosom/i and $skd='Mitochondrion';
		$skd='Mitochondrion';
	} elsif ( $stitle =~ m/\b(?:plastid)\b/i ) {
		$skd='Plastid';
	} elsif ( $stitle =~ m/\b(?:rDNA|rRNA|ribosomal RNA|ribosomal DNA)\b/i ) {
		$skd='rDNA';
	} elsif ( $skd eq 'Eukaryota' and $stitle =~ m/\b(?:Satellite)\b/i ) {
		$skd='Satellite';
	}
	my ($qs, $qe) = @ta[6,7]; # Query start and end
	my $qlen = $ta[12]; # Query length

	if ($skd =~ m/;/) {
		$skd =~ s/;?Eukaryota;?//; 
	}

	defined $skingdom{$skd} or die "[Err] Unknwn kingdom definition.\n[$ta[17]]\n$_\n";
	push( @{$qBlks{$qid}}, [$qs, $qe, $skd, $sid] ); # Any other information needed? $qs <= $qe ;
	defined $qLen{$qid} or $qLen{$qid} = $ta[12];
	for my $tp ($qs, $qe) {
		defined $dvd_site{$qid}{$tp} or $dvd_site{$qid}{$tp} = 1;
	}
}
&tsmsg("[Msg] file read in.\n");


# my %qUnit; # Or I will directly print them out to save memory usage.
print {$oFH} join("\t", qw/qseqid qlen qstart qend qspan KingdomCounts InExcludeCounts/)."\n";
if ( $is_joinInEx == 1 ) {
	print {$joinInExFH} join("\t", qw/qseqid qlen qstart qend qspan KingdomCounts InExcludeCounts/)."\n";
}


for my $qid (sort { $scf_order{$a} <=> $scf_order{$b} } keys %dvd_site) {
	&tsmsg("[Msg]Processing scaff [$qid]\n");
	## Make minimum block units for each query.
	my @qpoints = sort { $a<=>$b } keys %{ $dvd_site{$qid} }; # Sort block boundary sites.
	# $qpoints[0] > 1 and unshift(@qpoints, 1);
	# $qpoints[-1] < $qLen{$qid} and push(@qpoints, $qLen{$qid});
	my @qMinBlks ; # Initialize minimum blocks. [qstart, qend, {{skingdom=>nCount}}, {hitName=>1} ]
	if ( @qpoints == 0 ) {
		&stopErr("[Err] Failed to find dvd_site for [$qid]\n"); 
	} elsif ( @qpoints == 1 ) {
		# This only happens when qS == qE in blocks; 
		if ($qpoints[0] == 1) {
			push( @qMinBlks, [ 0.5, 1.5 ] ); 
		} else {
			push( @qMinBlks, [ 0.5, $qpoints[0]-0.5 ] ); 
		}
		push( @qMinBlks, [ $qpoints[0]-0.5, $qpoints[0]+0.5 ] ); 
	} else {
		# Initialize @qMinBlks
		if ($qpoints[0] > 1) {
			push(@qMinBlks, [0.5, $qpoints[0]-0.5]); # Unmapped region in 5'-end.
		}
		push(@qMinBlks, [$qpoints[0]-0.5, $qpoints[0]+0.5]); # Left-most point of mapping region.
		for (my $i=1; $i<=$#qpoints; $i++) {
			if ( $qpoints[$i-1]+1 < $qpoints[$i] ) {
				push(@qMinBlks, [$qpoints[$i-1]+0.5, $qpoints[$i]-0.5]);
			}
			push(@qMinBlks, [$qpoints[$i]-0.5, $qpoints[$i]+0.5]);
		}
	}
	$qMinBlks[-1][1] > $qLen{$qid}+0.5 and &stopErr("[Err] qMinBlks[-1][1] $qMinBlks[-1][1] > $qLen{$qid}+0.5\n");
	if ( $qMinBlks[-1][1] < $qLen{$qid}+0.5 ) {
		# There is some region left in the 3'-end. 
		push( @qMinBlks, [ $qMinBlks[-1][1], $qLen{$qid}+0.5 ] ); 
	}
	## Give type counts (contributions) for each unit in qMinBlks according to qBlks.
	my @qSrtBlks = sort { $a->[0]<=>$b->[0] || $a->[1]<=>$b->[1] || $skingdom{$a->[2]} <=> $skingdom{$b->[2]} } @{ $qBlks{$qid} }; # Sort query aligned blocks for following comparison.
	for (my $minI = 0; $minI < @qMinBlks; $minI++) {
		my $minE = $qMinBlks[$minI]; # [ qstart, qend, {{skingdom=>nCount}} ]
		BLOC:
		for (my $j=0; $j<@qSrtBlks; $j++) {
			my $srtE = $qSrtBlks[$j]; # [ qstart, qend, skingdom, hitName ]
			$srtE->[1] < $minE->[0] and next BLOC;
			$srtE->[0] > $minE->[1] and last BLOC;
			# Since blocks in qMinBlks are the smallest ones, any overlapping qBlks should cover all of this minE block.
			if ( $only1Hit ) {
				# If count sequence coverage by target sequence number instead of HSPs.
				defined $minE->[3]{ $srtE->[3] } or $minE->[2]{ $srtE->[2] } ++;
				$minE->[3]{ $srtE->[3] } ++;
			} else {
				$minE->[2]{ $srtE->[2] } ++;
			}
		}#End BLOC:for
	}#End

	# Save the final classification.
	# $qUnit{$qid} = [@qMinBlks];
	# Or directly print them out to save memory usage.
	my @combInEx;
	for my $tr1 (@qMinBlks) {
		# [qstart, qend, {{skingdom=>nCount}} ]
		my (@oskd, @oInEx);
		my (@j_oskd, @j_oInEx); # Edit here.
		if (!defined $tr1->[2]) {
			# Not mapped by any hit
			@oskd  = ("Un:0");
			@oInEx = ("In:0");
		}else{
			#
			# my @skds = sort { $tr1->[2]{$b} <=> $tr1->[2]{$a} || $skingdom{$a}<=>$skingdom{$b} } keys %{ $tr1->[2] };
			my @skds = sort { $skingdom{$a}<=>$skingdom{$b} } keys %{ $tr1->[2] };
			my %nInEx;
			for my $skd (@skds) {
				push( @oskd, join(":","$skd",$tr1->[2]{$skd}) );
				$nInEx{ $isInEx{$skd} } += $tr1->[2]{$skd};
			}
			for my $ie (sort { $nInEx{$b} <=> $nInEx{$a} } keys %nInEx) {
				push( @oInEx, join(":", $ie, $nInEx{$ie}) );
			}
		}
		# print $oFH join("\t", qw/qseqid qlen qstart qend qspan KingdomCounts InExcludeCounts/)."\n";
		print {$oFH} join("\t",
			$qid,
			$qLen{$qid},
			$tr1->[0]+0.5,
			$tr1->[1]-0.5,
			$tr1->[1]-$tr1->[0],
			join(";;", @oskd),
			join(";;", @oInEx)
		)."\n";
		if ($is_joinInEx == 1) {
			my ($kR1, $vR1, $kStr1) = &parse1(\@oskd);
			my ($kR2, $vR2, $kStr2) = &parse1(\@oInEx);
			if ( scalar(@combInEx) == 0 ) {
				# Initialize the first qMinBlk.
				push( @combInEx,
					[ 1,                   # qstart
					  $tr1->[1]-0.5,        # qend
					  [$kR1, $vR1, $kStr1], # [KingdomTypes_Ref, KingdomCounts_Ref, KingdomTypes_String]
					  [$kR2, $vR2, $kStr2]  # [InExcludeTypes_Ref, InExcludeCounts_Ref, InExcludeTypes_String]
					]
				);
			} elsif ( $kStr2 eq 'In' and $vR2->[0] == 0 ) {
				# the current block is unknown because it is unmapped. "In:0";
				push( @combInEx,
					[ $tr1->[0]+0.5,           # qstart
					  $tr1->[1]-0.5,           # qend
					  [$kR1, $vR1, $kStr1],    # [KingdomTypes_Ref, KingdomCounts_Ref, KingdomTypes_String]
					  [$kR2, $vR2, $kStr2]     # [InExcludeTypes_Ref, InExcludeCounts_Ref, InExcludeTypes_String]
					]
				);
			#} elsif ( scalar(@combInEx) == 1 and $combInEx[-1][3][2] eq 'In' and $combInEx[-1][3][1][0] == 0 and ($tr1->[0] - 0.5 <= $opts{maxUn}) ) {
			#	# This is the 2nd block and the 1st block is SHORT Unknown.
			#	# [4/21/2022] I don't want to merge the first unknown block into known.
			#	pop(@combInEx);
			#	push( @combInEx,
			#		[1,
			#		$tr1->[1]+0.5,
			#		[$kR1, $vR1, $kStr1],
			#		[$kR2, $vR2, $kStr2]
			#		]
			#	);
			} elsif ( $#combInEx > 0 and $combInEx[-1][3][2] eq 'In' and $combInEx[-1][3][1][0] == 0 and $combInEx[-1][1]-$combInEx[-1][0]+1 <= $opts{'maxUn'} and !($combInEx[-2][3][2] eq 'In' and $combInEx[-2][3][1][0] == 0) and $combInEx[-2][3][2] eq $kStr2 ) {
				# The previous block is short unmapped region and the one before the previous one is mapped.
				# Same Kingdom class, so renew the [-2] element and remove the previous unmapped short region.
				pop(@combInEx);
				&renewEle( $combInEx[-1] ,
					[ $tr1->[0]+0.5, $tr1->[1]-0.5, [$kR1, $vR1, $kStr1], [$kR2, $vR2, $kStr2] ]
				);
			} elsif ( !($combInEx[-1][3][2] eq 'In' and $combInEx[-1][3][1][0] == 0) and $combInEx[-1][3][2] eq $kStr2 ) {
				# The previous block and the current block are both aligned and share the same InEx class.
				&renewEle( $combInEx[-1] ,
					[ $tr1->[0]+0.5, $tr1->[1]-0.5, [$kR1, $vR1, $kStr1], [$kR2, $vR2, $kStr2] ]
				);
			} else {
				# Just push a new record.
				push( @combInEx,
					[ $tr1->[0]+0.5,
					  $tr1->[1]-0.5,
					  [$kR1, $vR1, $kStr1],
					  [$kR2, $vR2, $kStr2]
					]
				);
			}# End if ( scalar(@combInEx) == 0 )
		}# End if ($is_joinInEx == 1)
	}# for my $tr1 # Used to output query units.

	# For combine
	if ( $is_joinInEx == 1 ) {
		# I don't want to remove the unmapped region at the ends which are not clipped by two known blocks.
		for my $tr2 ( @combInEx ) {
			my (@sdk_str, @inEx_str);
			for (my $i=0; $i<@{$tr2->[2][0]}; $i++) {
				push( @sdk_str , join(":", $tr2->[2][0][$i], sprintf("%.2f", $tr2->[2][1][$i])) );
			}
			for (my $i=0; $i<@{$tr2->[3][0]}; $i++) {
				push( @inEx_str, join(":", $tr2->[3][0][$i], sprintf("%.2f", $tr2->[3][1][$i])) );
			}
			print {$joinInExFH} join("\t",
				$qid,
				$qLen{$qid},
				$tr2->[0],
				$tr2->[1],
				$tr2->[1] - $tr2->[0] + 1,
				join(";;", @sdk_str),
				join(";;", @inEx_str)
			)."\n";
		}# End for my $tr2
	}# End if ( $is_joinInEx == 1 )
}# for my $qid

if ( $is_joinInEx == 1 ) {
	close( $joinInExFH );
}

&tsmsg("[Rec]Over.\n");

#############################################################################
###############  End main             #######################################
#############################################################################

#############################################################################
###############  Subroutines          #######################################
#############################################################################

# Parse elements in @oskd and @oInEx
# Input : \@oskd / \@oInEx
# Output: (\@keys, \@values, "k1::k2::k3...")
sub parse1 {
	my $inR=shift;
	my (@tk, @tv);
	for (@$inR) {
		m/^([^:]+):([\d.]+)$/ or die "[Err]Wrong format for KingdomCounts/InExcludeCounts [$_]\n";
		push(@tk, $1);
		push(@tv, $2);
	}
	return(\@tk, \@tv, join(';;', @tk));
}# End sub parse1

sub joinVR {
	my ($e1, $e2, $size1, $size2, $sizeB) = @_; 
	my %k2v; 
	for (my $i=0; $i<@{$e1->[0]}; $i++) {
		$k2v{$e1->[0][$i]} = $e1->[1][$i] * $size1; 
	}
	for (my $i=0; $i<@{$e2->[0]}; $i++) {
		$k2v{$e2->[0][$i]} += ( $e2->[1][$i] * $size2 ); 
	}
	for my $k (keys %k2v) {
		$k2v{$k} = $k2v{$k}/$sizeB; 
	}
	my @k = sort { $skingdom{$a}<=>$skingdom{$b} } keys %k2v; 
	my @v = @k2v{@k}; 
	return( [\@k , \@v, join(';;', @k)] ); 
}# joinVR() 

# Renew the last element in array @combInEx
# input  : ($combInEx[-1], $toBeAdd_combInEx_ele)
# output : new $combInEx[-1] .
# In fact i do not think we need that return value, because the passed variable is an reference.
sub renewEle {
	my ($er1, $er2) = @_;
	my $size1 = $er1->[1] - $er1->[0] + 1;
	my $size2 = $er2->[1] - $er2->[0] + 1;
	my $new_s = ($er1->[0] > $er2->[0]) ? $er2->[0] : $er1->[0]; 
	my $new_e = ($er1->[1] < $er2->[1]) ? $er2->[1] : $er2->[1]; 

	$er1->[2] = &joinVR( $er1->[2], $er2->[2], $size1, $size2, $new_e-$new_s+1 ); 
	$er1->[3] = &joinVR( $er1->[3], $er2->[3], $size1, $size2, $new_e-$new_s+1 ); 
	$er1->[0] = $new_s; 
	$er1->[1] = $new_e; 

	return($er1);
}# End sub renewEle

#############################################################################
###############  Dropped subroutines          ###############################
#############################################################################

