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
use Getopt::Long; 
use fileSunhh; 

my %opts; 

GetOptions(\%opts, 
	"help!", 
	"protKLfile:s", "onlyComplete!", 
	"noAddTgt!", 
	"outFile:s", 
	"outAug!", 
	"scafKLfile:s", "trimOverflow!", 
); 

sub usage {
	print <<UU;
##########################################################################################
# perl $0 in.spaln.gff3
#
# -protKLfile     [filename] key\\tlength for protein. 
# -onlyComplete   Only output complete aligned models. 
#
# -noAddTgt       Not add target_name to mRNA element again. 
#
# -outAug         Output augustus format gff for training with gff2gbSmallDNA.pl. 
#
# -trimOverflow   Remove scaffold regions overflowing the maximum length of scaffold. 
# -scafKLfile     [filename] key\\tlength for scaffold. 
# 
# -outFile        [filename] file to create and write. 
##########################################################################################
UU
	exit 1; 
}

-t and !@ARGV and &usage(); 
$opts{'help'} and &usage(); 

my %k2l_scaf; 
if ( defined $opts{'scafKLfile'} or defined $opts{'trimOverflow'} ) {
	my $fh = &openFH( $opts{'scafKLfile'}, '<' ); 
	while (<$fh>) {
		chomp; m/^\s*$/ and next; 
		my @ta = split(/\t/, $_); 
		$k2l_scaf{$ta[0]} = $ta[1]; 
	}
	close($fh); 
}

my %k2l_prot; 
if ( defined $opts{'protKLfile'} or defined $opts{'onlyComplete'} ) {
	my $fh = &openFH( $opts{'protKLfile'}, '<' ); 
	while (<$fh>) {
		chomp; m/^\s*$/ and next; 
		my @ta = split(/\t/, $_); 
		$k2l_prot{$ta[0]} = $ta[1]; 
	}
	close($fh); 
}

my $oFH = ( defined $opts{'outFile'} ) ? &openFH($opts{'outFile'}, '>') : \*STDOUT ; 

my @geneLines; 
my $mID = ''; 
while (<>) {
	m/^#/ and next; 
	m/^\s*$/ and next; 
	m/^>/ and last; 
	chomp; 
	my @ta = split(/\t/, $_); 
	if ( $ta[2] eq 'mRNA' or $ta[2] eq 'protein_match' ) {
		if ( scalar(@geneLines) > 0 ) {
			$opts{'noAddTgt'} or &add_tgt(\@geneLines); 
			&rm_bad1st(\@geneLines); 
			$opts{'trimOverflow'} and &trim_overflow(\@geneLines, \%k2l_scaf, $opts{'trimOverflow'}); 
			&writeComplete($oFH, \@geneLines, \%k2l_prot, $opts{'onlyComplete'}); 
			@geneLines = (); 
			$mID = ''; 
		}
		if ( $opts{'outAug'} ) {
			$ta[8] =~ s!^ID=([^\s;]+)(?:;Parent=([^\s;]+))?(?:;Name=([^\s;]+))?;?$!transcript_id \"$1\"! or &stopErr("[Err] 2: $_\n"); 
			$mID = $1; 
			push(@geneLines, [[@ta], '']); 
		} else {
			$ta[8] =~ s!^ID=([^\s;]+)(?:;Parent=([^\s;]+))?(?:;Name=([^\s;]+))?;?$!ID=$1! or &stopErr("[Err] 2: $_\n"); 
			$ta[1] = "blastx"; 
			$ta[2] = "protein_match"; 
			push(@geneLines, [[@ta], '']); 
		}
	} elsif ( $ta[2] eq 'cds' or $ta[2] eq 'match_part' ) {
		if ( $opts{'outAug'} ) {
			$ta[8] =~ s!^(?:ID=[^\s;]+;)?Parent=([^\s;]+)(?:;Name=[^\s;]+)?(;Target=(\S+) (\d+) (\d+) [+-]);?$!transcript_id \"$1\"$2! or &stopErr("[Err] 5: $_\n"); 
			push(@geneLines, [[@ta], $3, [$4, $5] ]); 
		} else {
			$ta[8] =~ s!^(ID=[^\s;]+;Parent=[^\s;]+)(?:;Name=[^\s;]+)?(;Target=(\S+) (\d+) (\d+) [+-]);?$!$1$2! or &stopErr("[Err] 5: $_\n"); 
			$ta[1] = "blastx"; 
			$ta[2] = "match_part"; 
			push(@geneLines, [[@ta], $3, [$4, $5] ]); 
		}
	} else {
		$ta[2] eq 'gene' and next; 
		&stopErr("[Err] 1: $_\n"); 
	}
}

if ( scalar(@geneLines) > 0 ) {
	$opts{'noAddTgt'} or &add_tgt(\@geneLines); 
	&rm_bad1st(\@geneLines); 
	$opts{'trimOverflow'} and &trim_overflow(\@geneLines, \%k2l_scaf, $opts{'trimOverflow'}); 
	&writeComplete($oFH, \@geneLines, \%k2l_prot, $opts{'onlyComplete'});
	@geneLines = (); 
}

sub writeComplete {
	my ($fh, $ar, $hr_prot, $should_complete) = @_; 
	( defined $should_complete and $should_complete ) or $should_complete = 0; 

	my $tgt = $ar->[1][1]; 
	if ( $should_complete ) {
		( defined $hr_prot->{$tgt} and $hr_prot->{$tgt} > 0 ) or &stopErr("[Err] tgt [$tgt] length err.\n");
		my $is_complete = 1;

		my $prevE = 0;
		for ( my $i=1; $i<@$ar; $i++ ) {
			$ar->[$i][0][8] =~ m/;Target=$tgt (\d+) (\d+) [+\-][\s;]*$/ or &stopErr("[Err] 2: $ar->[$i][0][8]\n");
			my ($curS, $curE) = ($1, $2);
			$curS == $prevE + 1 or do { $is_complete = 0; last; } ;
			$prevE = $curE;
		}
		$prevE == $hr_prot->{$tgt} or $is_complete = 0;
		$is_complete == 0 and return; 
	}

	for (my $i=0; $i<@$ar; $i++) {
		print {$fh} join("\t", @{$ar->[$i][0]})."\n";
	}
	return 0;
}

sub trim_overflow {
	my ($ar, $hr_scaf, $should_trim) = @_; 
	( defined $should_trim and $should_trim ) or $should_trim = 0; 

	$should_trim == 0 and return 0;

	my ($gen_s, $gen_e); 
	my $scaf_id = $ar->[0][0][0]; 
	defined $hr_scaf->{$scaf_id} or &stopErr("[Err] No length for scaffold [$scaf_id]\n"); 
	my $scaf_len = $hr_scaf->{$scaf_id}; 
	for (my $i=1; $i<@$ar; $i++) {
		if ( $ar->[$i][0][3] > $scaf_len ) {
			&tsmsg("[Err] We met a bad element!\n"); 
			&tsmsg("[Err] scaff_len{$scaf_id} = $scaf_len\n"); 
			&tsmsg("[Err] ", join("\t", @{$ar->[$i][0]}), "\n"); 
			&stopErr("[Err] \n"); 
		} elsif ( $ar->[$i][0][4] > $scaf_len ) {
			$ar->[$i][0][4] = $scaf_len ; 
		}
		defined $gen_s or $gen_s = $ar->[$i][0][3]; 
		defined $gen_e or $gen_e = $ar->[$i][0][4]; 
		$gen_s < $ar->[$i][0][3] and $gen_s = $ar->[$i][0][3]; 
		$gen_e > $ar->[$i][0][4] and $gen_e = $ar->[$i][0][4]; 
	}
	return 0; 
}# trim_overflow

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

