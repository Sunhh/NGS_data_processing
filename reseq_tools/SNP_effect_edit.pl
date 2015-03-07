#!/usr/bin/perl -w 
# 方法一: 20120110 不考虑邻近SNP之间的互相影响，独立检查SNP所带来的变化;
#         SNP影响分类: 
#            1. 基因间区 : 1.1 基因上游500bp 1.2 基因下游500bp 1.3 其余基因间区; 
#            2. 基因内含子区: 2.1 intron/exon边界(即exon边界碱基外延2bp) 2.2 其它内含子区域; 
#            3. Coding区: 3.1 起始密码子变为非起始密码子  3.2 终止密码子变化为继续编码  
#                         3.3 氨基酸编码变化-转义  3.4 氨基酸编码变化-终止   3.5 同义突变; 
#                         3.6 终止密码子替换为其它终止密码子; 3.7 起始密码子变为其它起始密码(AA变化); 
# edit20120118:           3.8 位于编码区, 但氨基酸编码中存在'N', 因而无法判断; 
# ClV6 genome intron length stat: SUM             MEAN                    MEDIAN  MIN     MAX     Count   NoNull
#                                 38918967        464.338157392383        199     11      9976    83816   83816
# 2013-08-19 开始遇到太多Indel、杂合等问题，这些单独处理太麻烦，而且有时基因位置信息也是有用的。

use strict; 
-t and !@ARGV and die "perl $0 in_SNP\n"; 

my %codon; 
# Initial %codon; 
{
	my @aa = split(//, 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'); 
	my @ss = split(//, '---M---------------M---------------M----------------------------'); 
	my @l1 = split(//, 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'); 
	my @l2 = split(//, 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'); 
	my @l3 = split(//, 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'); 
	for (my $i=0; $i<@aa; $i++) {
		my $setag = '-'; 
		$aa[$i] eq '*' and $setag = '*'; 
		$ss[$i] eq 'M' and $setag = 'M'; 
		my $bbb = join('', $l1[$i], $l2[$i], $l3[$i]); 
		$codon{$bbb} = [$setag, $aa[$i]]; 
	}
}

#先偷个懒, 直接读取CDS序列; 
my $cdsfile = '/share/app/watermelon/share/V6_WM_20100930/GeneAnnotation/V6_WM_CDS'; 
my $gfffile = '/share/app/watermelon/share/V6_WM_20100930/GeneAnnotation/watermelon_v6.scaffold.glean.1015.gff.150_filter.gff_CHRloc'; 


$cdsfile = '/data/Sunhh/database/Watermelon/WM97/WM97_v6.cds.fasta';
$gfffile = '/data/Sunhh/database/Watermelon/WM97/WM97_v6.annot.gff3';

open CDS,'<',"$cdsfile" or die; 
my (%cdsseq, $tkey); 
while (<CDS>) {
	s/\s+$//; 
	/^\s*$/ and next; 
	if (/^>(\S+)/) {
		$tkey = $1; 
	}else{
		s/\s//g; 
		$cdsseq{$tkey} .= $_; 
	}
}
close CDS; 

open GFF,'<',"$gfffile" or die; 
my %anno; 
while (<GFF>) {
	chomp; s/\s+$//; 
	/^\s*$/ and next; 
	my @ta = split(/\t/, $_); 
	my ($cid, $t1, $ts, $te, $tstrand, $tname) = @ta[0, 2, 3, 4, 6, 8]; 
	if ($t1 eq 'mRNA') {
		$tname =~ s/^ID=([^;=\s]+);$/$1/ or die "Failed $tname\n"; 
		defined $cdsseq{$tname} or die "No cdsseq for $tname\n"; 
		push(@{$anno{$cid}}, [[$ts,$te], $tstrand, [], $tname, length($cdsseq{$tname}), $cdsseq{$tname}]); 
	}elsif ($t1 eq 'CDS') {
		$tname =~ s/^Parent=([^;=\s]+);$/$1/ or die "Failed1 $tname\n"; 
		$anno{$cid}[-1][3] eq $tname or die "Diff GenID in: $anno{$cid}[-1][3] eq $tname\n"; 
		push(@{$anno{$cid}[-1][2]}, [$ts, $te]); 
	}else{
		$t1 eq 'gene' and next; 
		die "Line:$_\n$t1\n"; 
	}
}
close GFF; 

SNP_LINE: 
while (<>) {
	chomp; s/\s+$//; 
	/^\s*$/ and next; 
#	/^chr(?:omosome|om)?\t/i and next; 
	my @ta = split(/\t/, $_); 
	my ($chr_id, $chr_pos) = @ta[0,1]; 
#	my $refbase = $ta[2]; $refbase = uc($refbase); 
	my $refbase = undef(); 
	my $newbase = undef(); 
	
	
	my @td; # Storing the effections. 
	if ($ta[0] =~ m/^chr(?:omosome|om)?$/) {
		print STDOUT join("\t", "SNP_Effect", "$_")."\n"; 
		next SNP_LINE; 
	} else {
		$refbase = $ta[2]; $refbase = uc($refbase); 
		for my $tb (@ta[3..$#ta]) {
			$tb = uc($tb); 
			$tb eq 'N' and next; 
#			length($tb) > 1 and do { warn "skip base [$tb]\n"; next; }; 
			$tb ne $refbase and do {$newbase = $tb; last; }; 
		}
		if ($newbase eq '*') {
			push(@td, 'Deletion'); 
			$newbase = $refbase; 
		} elsif ($newbase =~ m/\+/) {
			push(@td, 'Insertion'); 
			$newbase = $refbase; 
		} elsif ($newbase =~ m/\*/ ) {
			push(@td, 'WithDeletion'); 
			$newbase = $refbase; 
		} elsif ($newbase =~ m/^[ATGC]{2,}$/) {
			push(@td, 'Heterozygous'); 
			$newbase = $refbase; 
		} elsif ($newbase =~ m/^[ATGC]$/) {
#			my @tc = &parseSNP('', $anno{$chr_id}, $chr_pos, $newbase, $refbase); 
#			for my $tcr (@tc) {
#				push(@td, join(":", @$tcr)); 
#			}
		} else {
			warn "[Warn]Failed to get allele:$_\n"; 
			$newbase = $refbase; 
			@td = ("0.0:Failed"); 
		}
	}#End if ($ta[0]) 
	
	my @tc = &parseSNP('', $anno{$chr_id}, $chr_pos, $newbase, $refbase); 
	for my $tcr (@tc) {
		push(@td, join(":", @$tcr)); 
	}
	
	print STDOUT join("\t", join(',', @td), $_)."\n"; 
}


# my ($cdsseqR, $gffR, $position, $newbase, $refbase) = @_; 
# gffR : [mRNA_S_E, Strand, CDS_SEs, GenID, cds_len, cds_seq]
# back : [effect_type1, effect_type2, ...]; 
# 显然, newbase是'+'链上的; 
sub parseSNP {
	my ($cdsseqR, $gffR, $position, $newbase, $refbase) = @_; 
	my $up_len = 500; 
	$up_len = 2000; # for 10kb upstream of genes; 
	my $down_len = 500; 
	$down_len = 2000; # for 10kb downstream of genes; 
	my $span = 2; 
	my @type; 
	GENE: 
	for my $r1 (@$gffR) {
		if ($r1->[1] eq '+') {
			$position < $r1->[0][0]-$up_len and next GENE; 
			$position > $r1->[0][1]+$down_len and next GENE; 
			if ($position < $r1->[0][0]) {
				push(@type, [$r1->[3], '1.1']); next GENE; 
			}elsif ($position > $r1->[0][1]) {
				push(@type, [$r1->[3], '1.2']); next GENE; 
			}else{
				# 在起止exon限定区间内; 
				my $cds_pos = 0; 
				for my $r2 (@{$r1->[2]}) {
					if ($position > $r2->[1]) {
						# intron2; 
						if ($position <= $r2->[1]+$span) {
							push(@type, [$r1->[3], '2.1']); next GENE; 
						}else{
							# else we should check next exon. 
							$cds_pos += ($r2->[1]-$r2->[0]+1); 
						}
					}elsif ($position < $r2->[0]) {
						# intron1; 
						if ($position >= $r2->[0]-$span) {
							push(@type, [$r1->[3], '2.1']); next GENE; 
						}else{ 
							# 此处需要配合 $position > $r2->[1] 后的 $position <= $r2->[1]+$span 共同使用才合适! 
							push(@type, [$r1->[3], '2.2']); next GENE; 
						}
					}else{
						# position fall in this exon now!
						$cds_pos += ($position-$r2->[0]+1); 
						# 判断cds_pos的类型; 
						my $tbase = $newbase; 
						my $ref_tbase = $refbase; 
						my @ttype = &chk_cds_pos($cds_pos, $tbase, \$r1->[5], $r1->[4], $ref_tbase); 
						push(@type, [$r1->[3], @ttype]); next GENE; 
					}
				}# end for 
				scalar(@type) > 0 or die "Error here!\n@$r1\n"; 
			}
		}elsif ($r1->[1] eq '-') {
			$position < $r1->[0][0]-$down_len and next GENE; 
			$position > $r1->[0][1]+$up_len and next GENE; 
			if ($position < $r1->[0][0]) {
				push(@type, [$r1->[3], '1.2']); next GENE; 
			}elsif ($position > $r1->[0][1]) {
				push(@type, [$r1->[3], '1.1']); next GENE; 
			}else{
				# 在起止exon限定区间内; 
				my $cds_pos = 0; 
				for my $r3 (@{$r1->[2]}) {
					if ($position < $r3->[0]) {
						# intron2; 
						if ($position >= $r3->[0]-$span) {
							push(@type, [$r1->[3], '2.1']); next GENE; 
						}else{
							$cds_pos += ($r3->[1]-$r3->[0]+1); 
						}
					}elsif ($position > $r3->[1]) {
						# intron1; 
						if ($position <= $r3->[1]+$span) {
							push(@type, [$r1->[3], '2.1']); next GENE; 
						}else{
							push(@type, [$r1->[3], '2.2']); next GENE; 
						}
					}else{
						# position fall in this exon now!
						$cds_pos += ($r3->[1]-$position+1); 
						my $tbase = $newbase; 
						$tbase =~ tr/ATGC/TACG/; 
						my $ref_tbase = $refbase; 
						$ref_tbase =~ tr/ATGC/TACG/; 
						my @ttype = &chk_cds_pos($cds_pos, $tbase, \$r1->[5], $r1->[4], $ref_tbase); 
						push(@type, [$r1->[3], @ttype]); next GENE; 
					}
				}# end for 
				scalar(@type) > 0 or die "Err here!\n"; 
			}
		}else{ # Failed to know '+'/'-'
			die "Failed to parse strand for $r1->[1]!\n"; 
		}
	}# end for gffR
	scalar(@type) == 0 and push(@type, ['NULL', '1.3']); 
	return (@type); 
}# end sub parseSNP 

# gffR : [mRNA_S_E, Strand, CDS_SEs, GenID, cds_len, cds_seq]
# back : [effect_type1, effect_type2, ...]; 
# 方法一: 20120110 不考虑邻近SNP之间的互相影响，独立检查SNP所带来的变化;
#         SNP影响分类: 
#            1. 基因间区 : 1.1 基因上游500bp 1.2 基因下游500bp 1.3 其余基因间区; 
#            2. 基因内含子区: 2.1 intron/exon边界(即exon边界碱基外延2bp) 2.2 其它内含子区域; 
#            3. Coding区: 3.1 起始密码子变为非起始密码子  3.2 终止密码子变化为继续编码  
#                         3.3 氨基酸编码变化-转义  3.4 氨基酸编码变化-终止   3.5 同义突变; 
#                         3.6 终止密码子替换为其它终止密码子; 3.7 起始密码子变为其它起始密码(AA变化); 
# ClV6 genome intron length stat: SUM             MEAN                    MEDIAN  MIN     MAX     Count   NoNull
#                                 38918967        464.338157392383        199     11      9976    83816   83816

# mut应该是与seq同链的碱基类型; 
sub chk_cds_pos{
	my ($p, $mut, $seqR, $len, $ref) = @_; 
	my ($newb3, $rawb3); 
	if ($p <= 3) {
		my $aa_pos = 1; 
		$rawb3 = $newb3 = substr($$seqR, 0, 3); 
		substr($newb3, $p-1, 1) = $mut; 
		substr($rawb3, $p-1, 1) = $ref; 
		if ($codon{$newb3}[0] eq '*') {
			return ('3.4', "$codon{$rawb3}[1]$aa_pos$codon{$newb3}[1]"); 
		}elsif ($codon{$newb3}[0] eq '-') {
			return ('3.1', "$codon{$rawb3}[1]$aa_pos$codon{$newb3}[1]"); 
		}elsif ($codon{$newb3}[1] ne $codon{$rawb3}[1]) {
			return ('3.7', "$codon{$rawb3}[1]$aa_pos$codon{$newb3}[1]"); 
		}else{
			return ('3.5', "$codon{$rawb3}[1]$aa_pos$codon{$newb3}[1]"); 
		}
#	}elsif ($p > $len-3) {
#		my $aa_pos = int($len/3); 
#		$newb3 = $rawb3 = substr($$seqR, $len-3, 3); 
#		$codon{$rawb3}[0] eq '*' or die "Tail BBB is $rawb3, not * but $codon{$rawb3}[1]\n"; 
#		substr($newb3, 2-($len-$p), 1) = $mut; 
#		substr($rawb3, 2-($len-$p), 1) = $ref; 
#		if ($codon{$newb3}[0] ne '*') {
#			return ('3.2', "*$aa_pos$codon{$newb3}[1]"); 
#		}else{
#			return ('3.6'); 
#		}
	}else{
		my $aa_pos = int($p/3)+1; 
		my $frame = $p%3; 
		if ($frame == 0) {
			$frame = 3; 
			$aa_pos--; 
		}
		my $ssp = $p-$frame+1; 
		$rawb3 = $newb3 = substr($$seqR, $ssp-1, 3); 
		(defined $newb3 and $newb3 ne '') or die "$ssp-1, $$seqR\n"; 
		substr($newb3, $frame-1, 1) = $mut; 
		substr($rawb3, $frame-1, 1) = $ref; 

		if (!(defined $codon{$newb3}[0] and defined $codon{$rawb3})) {
			if ($newb3 =~ /N/i or $rawb3 =~ /N/i) {
				my $rawX = (defined $codon{$rawb3}) ? $codon{$rawb3} : 'X' ; 
				my $newX = (defined $codon{$newb3}) ? $codon{$newb3} : 'X' ; 
				return ('3.8', "$rawX${aa_pos}$newX"); 
			}else{
				warn "Pos=$p, Mut=$mut, RawBBB=$rawb3, NewBBB=$newb3\n"; 
			}
		}
		if ($codon{$rawb3}[0] ne '*') { 
			if ($codon{$newb3}[0] eq '*') {
				return ('3.4', "$codon{$rawb3}[1]$aa_pos$codon{$newb3}[1]"); 
			}elsif ($codon{$newb3}[1] ne $codon{$rawb3}[1]) {
				return ('3.3', "$codon{$rawb3}[1]$aa_pos$codon{$newb3}[1]"); 
			}else{
				return ('3.5', "$codon{$rawb3}[1]$aa_pos$codon{$newb3}[1]"); 
			}
		}else{
			if ($codon{$newb3}[0] ne '*') {
				return ('3.2', "*$aa_pos$codon{$newb3}[1]"); 
			} else {
				return ('3.6'); 
			}
		}#End if $codon{$rawb3}[0] ne '*'
	}
}# end sub chk_cds_pos
