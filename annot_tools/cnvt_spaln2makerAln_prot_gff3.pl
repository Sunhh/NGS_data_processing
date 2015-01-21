#!/usr/bin/perl
# 20150119 : I found spaln somtimes gave bad exon result: 
# ##sequence-region       S401681_pilon 58369 90112
# S401681_pilon   ALN     gene    58573   72643   859     +       .       ID=gene00001;Name=S401681_pilon_65
# S401681_pilon   ALN     mRNA    58573   72643   859     +       .       ID=mRNA00001;Parent=gene00001;Name=S401681_pilon_65
# S401681_pilon   ALN     cds     58573   58580   -150    +       0       ID=cds00001;Parent=mRNA00001;Name=S401681_pilon_65;Target=Bv_16390_kjfx.t1 1 0 +
# S401681_pilon   ALN     cds     71401   72643   1085    +       1       ID=cds00002;Parent=mRNA00001;Name=S401681_pilon_65;Target=Bv_16390_kjfx.t1 1 395 +
# S402351_pilon   ALN     gene    8409    25030   1645    +       .       ID=gene14902;Name=S402351_pilon_16
# S402351_pilon   ALN     mRNA    8409    25030   1645    +       .       ID=mRNA14902;Parent=gene14902;Name=S402351_pilon_16
# S402351_pilon   ALN     cds     8409    8412    -138    +       0       ID=cds83202;Parent=mRNA14902;Name=S402351_pilon_16;Target=Bv8_194420_dcea.t1 1 0 +
# S402351_pilon   ALN     cds     20587   20679   110     +       2       ID=cds83203;Parent=mRNA14902;Name=S402351_pilon_16;Target=Bv8_194420_dcea.t1 1 26 +
# S402351_pilon   ALN     cds     20768   20871   166     +       2       ID=cds83204;Parent=mRNA14902;Name=S402351_pilon_16;Target=Bv8_194420_dcea.t1 27 61 +
# S402351_pilon   ALN     cds     21303   21502   290     +       0       ID=cds83205;Parent=mRNA14902;Name=S402351_pilon_16;Target=Bv8_194420_dcea.t1 62 128 +
# S402351_pilon   ALN     cds     23070   23170   163     +       1       ID=cds83206;Parent=mRNA14902;Name=S402351_pilon_16;Target=Bv8_194420_dcea.t1 129 161 +
# S402351_pilon   ALN     cds     23274   23466   284     +       2       ID=cds83207;Parent=mRNA14902;Name=S402351_pilon_16;Target=Bv8_194420_dcea.t1 162 226 +
# S402351_pilon   ALN     cds     23559   23720   302     +       1       ID=cds83208;Parent=mRNA14902;Name=S402351_pilon_16;Target=Bv8_194420_dcea.t1 227 280 +
# S402351_pilon   ALN     cds     23828   23909   140     +       1       ID=cds83209;Parent=mRNA14902;Name=S402351_pilon_16;Target=Bv8_194420_dcea.t1 281 307 +
# S402351_pilon   ALN     cds     24290   24556   397     +       0       ID=cds83210;Parent=mRNA14902;Name=S402351_pilon_16;Target=Bv8_194420_dcea.t1 308 396 +
# S402351_pilon   ALN     cds     24653   24780   191     +       0       ID=cds83211;Parent=mRNA14902;Name=S402351_pilon_16;Target=Bv8_194420_dcea.t1 397 440 +
# S402351_pilon   ALN     cds     24883   25030   273     +       1       ID=cds83212;Parent=mRNA14902;Name=S402351_pilon_16;Target=Bv8_194420_dcea.t1 441 489 +
#
use strict; 
use warnings; 
use LogInforSunhh; 

-t and !@ARGV and die "perl $0 in.spaln.gff3\n"; 

my @geneLines; 
while (<>) {
	m/^#/ and next; 
	m/^\s*$/ and next; 
	m/^>/ and last; 
	chomp; 
	my @ta = split(/\t/, $_); 
	if ( $ta[2] eq 'mRNA' ) {
		if ( scalar(@geneLines) > 0 ) {
			&add_tgt(\@geneLines); 
			&rm_bad1st(\@geneLines); 
			for my $ar (@geneLines) {
				print STDOUT join("\t", @{$ar->[0]}) . "\n"; 
			}
			@geneLines = (); 
		}
		$ta[8] =~ s!^ID=([^\s;]+);Parent=([^\s;]+);Name=([^\s;]+)$!ID=$1! or &stopErr("[Err] 2: $_\n"); 
		$ta[1] = "blastx"; 
		$ta[2] = "protein_match"; 
		push(@geneLines, [[@ta], '']); 
	} elsif ( $ta[2] eq 'cds' ) {
		$ta[8] =~ s!^(ID=[^\s;]+;Parent=[^\s;]+;)Name=[^\s;]+;(Target=(\S+) (\d+) (\d+) [+-])$!$1$2! or &stopErr("[Err] 5: $_\n"); 
		$ta[1] = "blastx"; 
		$ta[2] = "match_part"; 
		push(@geneLines, [[@ta], $3, [$4, $5] ]); 
	} else {
		$ta[2] eq 'gene' and next; 
		&stopErr("[Err] 1: $_\n"); 
	}
}

if ( scalar(@geneLines) > 0 ) {
	&add_tgt(\@geneLines); 
	&rm_bad1st(\@geneLines); 
	for my $ar (@geneLines) {
		print STDOUT join("\t", @{$ar->[0]}) . "\n"; 
	}
	@geneLines = (); 
}

sub add_tgt {
	my $ar = shift; 
	my $tgt = ''; 
	for (my $i=1; $i<@$ar; $i++) {
		$tgt eq '' and $tgt = $ar->[$i][1]; 
		$tgt ne '' or &stopErr("[Err] 3: bad tgt [$tgt] in $ar->[0]\n"); 
		$tgt eq $ar->[$i][1] or &stopErr("[Err] 4: diff tgt [$tgt] vs. $ar->[$i][1] in $ar->[0]\n"); 
	}
	$ar->[0][0][8] .= ";Name=$tgt"; 
	return 0; 
}
sub rm_bad1st {
	my $ar = shift; 
	if ( $ar->[1][2][1] == 0 or $ar->[1][2][1] < $ar->[1][2][0] ) {
		my ($min, $max) = (-1, -1); 
		for ( my $i=2; $i<@$ar; $i++ ) {
			$min == -1 and $min = $ar->[$i][0][3]; 
			$max == -1 and $max = $ar->[$i][0][4]; 
			$min > $ar->[$i][0][3] and $min = $ar->[$i][0][3]; 
			$max < $ar->[$i][0][4] and $max = $ar->[$i][0][4]; 
		}
		$ar->[0][0][3] = $min; 
		$ar->[0][0][4] = $max; 
		&tsmsg("[Wrn] Remove alignment element: ", join("\t", @{$ar->[1][0]}), "\n"); 
		splice(@$ar, 1, 1); 
	}
	scalar(@$ar) >= 2 or &stopErr("[Err] No good alignment for : ", join("\t", @{$ar->[0][0]}), "\n"); 
	return 0; 
}

