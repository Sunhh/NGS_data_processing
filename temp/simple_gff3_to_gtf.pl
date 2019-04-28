#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

-t and !@ARGV and die "perl $0 WM97_v1.chr.annot_prot.gff3 > WM97_v1.chr.annot_prot.gtf\n"; 

my %l1_gene; 
my %gene_to_mRNA; 
my %l2_mRNA; 
my %l3_CDS; 
my %l3_exon; 

while (<>) {
	m!^\s*#|^\s*$! and next; 
	chomp; 
	my @ta=split(/\t/, $_); 
	# Split $ta[8]; 
	my %attr; 
	my $tmp_ta8 = $ta[8]; 
	while ($tmp_ta8 =~ s!^\s*([^=\s]+)\s*=\"?([^;"]+)\"?(?:;|\s*$)!!) {
		my ($k, $v) = ($1, $2); 
		$v =~ s!^\s+|\s+$!!g; 
		$attr{uc($k)} = $v; 
	}
	$tmp_ta8 =~ m!^\s*$! or die "Bad rest ta8 [$tmp_ta8] in line: $_\n"; 
	if      ( $ta[2] =~ m!^gene$!i ) {
		defined $attr{'ID'} or die "ID?\n$_\n"; 
		$attr{'NAME'} //= $attr{'ID'}; 
		my $new_ta8 = "gene_id \"$attr{'ID'}\"; gene_name \"$attr{'NAME'}\""; 
		for my $k1 (sort keys %attr) {
			$k1 =~ m!^(ID|NAME)$!i and next; 
			$new_ta8 .= " $k1 \"$attr{$k1}\";"; 
		}
		# print STDOUT join("\t", $ta[0], $ta[1], "gene", @ta[3..7], $new_ta8)."\n"; 
		defined $l1_gene{ $attr{'ID'} } and die "repG: $_\n"; 
		$l1_gene{$attr{'ID'}} = [ $ta[0], $ta[1], "gene", @ta[3..7], $new_ta8 ]; 
	} elsif ( $ta[2] =~ m!^(mRNA|transcript)$!i ) {
		defined $attr{'ID'} or die "ID?\n$_\n"; 
		defined $attr{'PARENT'} or die "Parent?\n$_\n"; 
		my $new_ta8 = "transcript_id \"$attr{'ID'}\"; gene_id \"$attr{'PARENT'}\";"; 
		for my $k1 (sort keys %attr) {
			$k1 =~ m!^(ID|PARENT)$!i and next; 
			$new_ta8 .= " $k1 \"$attr{$k1}\";"; 
		}
		defined $l2_mRNA{$attr{'ID'}} and die "repM: $_\n"; 
		$l2_mRNA{$attr{'ID'}} = [ $ta[0], $ta[1], "transcript", @ta[3..7], $new_ta8 ]; 
		push(@{$gene_to_mRNA{$attr{'PARENT'}}}, $attr{'ID'}); 
	} elsif ( $ta[2] =~ m!^exon$!i ) {
		defined $attr{'PARENT'} or die "Parent?\n$_\n"; 
		my $new_ta8 = "transcript_id \"$attr{'PARENT'}\";"; 
		for my $k1 (sort keys %attr) {
			$k1 =~ m!^(PARENT)$!i and next; 
			if ($k1 eq 'ID') {
				$new_ta8 .= " exon_id \"$attr{'ID'}\";"; 
			} else {
				$new_ta8 .= " $k1 \"$attr{$k1}\";"; 
			}
		}
		push(@{$l3_exon{$attr{'PARENT'}}}, [ $ta[0], $ta[1], "exon", @ta[3..7], $new_ta8 ]); 
	} elsif ( $ta[2] =~ m!^cds$!i ) {
		defined $attr{'PARENT'} or die "Parent?\n$_\n"; 
		my $new_ta8 = "transcript_id \"$attr{'PARENT'}\";"; 
		for my $k1 (sort keys %attr) {
			$k1 =~ m!^(PARENT)$!i and next; 
			if ($k1 eq 'ID') {
				$new_ta8 .= " cds_id \"$attr{'ID'}\";"; 
			} else {
				$new_ta8 .= " $k1 \"$attr{$k1}\";"; 
			}
		}
		push(@{$l3_CDS{$attr{'PARENT'}}}, [ $ta[0], $ta[1], "cds", @ta[3..7], $new_ta8 ]); 
	} else {
		&tsmsg("[Wrn] Skip line for [$ta[2]] : $_\n"); 
	}
}

for my $gid (sort {$l1_gene{$a}[0] cmp $l1_gene{$b}[0] || $l1_gene{$a}[3] <=> $l1_gene{$b}[3] || $l1_gene{$a}[4] <=> $l1_gene{$a}[4]} keys %l1_gene) {
	print STDOUT join("\t", @{$l1_gene{$gid}})."\n"; 
	# print STDOUT join("\t", @{$l1_gene{$gid}}[0..7], "gene_id \"$gid\"; gene_type \"protein_coding\";")."\n"; 
	defined $gene_to_mRNA{$gid} or next; 
	my $genename = $gid; 
	$l1_gene{$gid}[8] =~ m!gene_name \"([^"]+)\"! and $genename = $1; 
	for my $mid (sort { $l2_mRNA{$a}[3] <=> $l2_mRNA{$b}[3] || $l2_mRNA{$a}[4] <=> $l2_mRNA{$b}[4] } @{$gene_to_mRNA{$gid}}) {
		$l2_mRNA{$mid}[8] =~ m/gene_id \"/ or $l2_mRNA{$mid}[8] .= " gene_id \"$gid\";"; 
		$l2_mRNA{$mid}[8] =~ m/gene_name \"/ or $l2_mRNA{$mid}[8] .= " gene_name \"$genename\";"; 
		print STDOUT join("\t", @{$l2_mRNA{$mid}})."\n"; 
		# $l2_mRNA{$mid}[2] = 'transcript'; 
		# print STDOUT join("\t", @{$l2_mRNA{$mid}}[0..7], "transcript_id \"$mid\"; transcript_type \"protein_coding\"; gene_id \"$gid\"; gene_type \"protein_coding\";")."\n"; 
		if ( defined $l3_exon{$mid} ) {
			if ($l2_mRNA{$mid}[6] eq '+') {
				@{$l3_exon{$mid}} = sort { $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4] } @{$l3_exon{$mid}}; 
			} elsif ($l2_mRNA{$mid}[6] eq '-') {
				@{$l3_exon{$mid}} = sort { $b->[4] <=> $a->[4] || $b->[3] <=> $a->[4] } @{$l3_exon{$mid}}; 
			} else {
				die "@{$l2_mRNA{$mid}}\n"; 
			}
			my $en = 0; 
			for my $cline (@{$l3_exon{$mid}}) {
				$cline->[8] =~ m/gene_id \"/   or $cline->[8] .= " gene_id \"$gid\";"; 
				$cline->[8] =~ m/gene_name \"/ or $cline->[8] .= " gene_name \"$genename\";"; 
				print STDOUT join("\t", @{$cline})."\n"; 
				$en ++; 
				# print STDOUT join("\t", @{$cline}[0..7], "exon_number $en; transcript_id \"$mid\"; transcript_type \"protein_coding\"; gene_id \"$gid\"; gene_type \"protein_coding\";")."\n"; 
			}
		}
		if ( defined $l3_CDS{$mid} ) {
			if ($l2_mRNA{$mid}[6] eq '+') {
				@{$l3_CDS{$mid}} = sort { $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4] } @{$l3_CDS{$mid}}; 
			} elsif ($l2_mRNA{$mid}[6] eq '-') {
				@{$l3_CDS{$mid}} = sort { $b->[4] <=> $a->[4] || $b->[3] <=> $a->[4] } @{$l3_CDS{$mid}}; 
			} else {
				die "@{$l2_mRNA{$mid}}\n"; 
			}
			my $en = 0; 
			for my $cline (@{$l3_CDS{$mid}}) {
				$cline->[8] =~ m/gene_id \"/   or $cline->[8] .= " gene_id \"$gid\";"; 
				$cline->[8] =~ m/gene_name \"/ or $cline->[8] .= " gene_name \"$genename\";"; 
				print STDOUT join("\t", @{$cline})."\n"; 
				$en ++; 
				# print STDOUT join("\t", @{$cline}[0..7], "exon_number $en; transcript_id \"$mid\"; transcript_type \"protein_coding\"; gene_id \"$gid\"; gene_type \"protein_coding\";")."\n"; 
			}
		}
	}
}


