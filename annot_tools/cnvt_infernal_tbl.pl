#!/usr/bin/perl
# 2019-05-07 Fix for target of description. 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"outBed!", 
	"outGff!", 
	"noTrunc!", 
); 


my $usage = <<HH; 

perl $0 P1Genom_V1.scf.fa.rRNA.tbl > P1Genom_V1.scf.fa.rRNA.tbl.tab

-outBed
-outGff

-noTrunc

HH

-t and !@ARGV and &LogInforSunhh::usage($usage); 
$opts{'help'} and &LogInforSunhh::usage($usage);

my @aln_tbl = &load_cmscan_tbl_1(); 
if ( $opts{'outBed'} ) {
	&out_file(\@aln_tbl, 'outBed'); 
} elsif ( $opts{'outGff'} ) {
	&out_file(\@aln_tbl, 'outGff'); 
} else {
	&out_file(\@aln_tbl, 'outTbl'); 
}

sub out_file {
	my ($lines_aref, $outType) = @_; 
	$outType //= 'outTbl'; 
	
	my %coln; 
	# $coln{'seq_from'} = 7; 
	# $coln{'seq_to'}   = 8; 
	# $coln{'score'}    = 14; 
	# $coln{'strand'}   = 9; 
	# $coln{'bias'}     = 13; 
	# $coln{'target_name'}  = 0; 
	# $coln{'E-value'}  = 15; 
	# $coln{'inc'}      = 16; 
	# $coln{'pass'}     = 11; 
	# $coln{'trunc'}    = 10; 
	# $coln{'query_name'}= 2; 
	my $hd = shift(@$lines_aref); 
	for (my $j=0; $j<@$hd; $j++) {
		if (defined $coln{$hd->[$j]}) {
			if ($hd->[$j] eq 'accession') {
				if      ($hd->[$j-1] eq 'target_name') {
					$hd->[$j] = 'target_accession'; 
				} elsif ($hd->[$j-1] eq 'query_name') {
					$hd->[$j] = 'query_accession'; 
				} else {
					die "repeat header 1 $hd->[$j]\n"; 
				}
			} else {
				die "repeat header 2 $hd->[$j]\n"; 
			}
		}
		defined $coln{$hd->[$j]} and die "repeat header 3 $hd->[$j]\n"; 
		$coln{$hd->[$j]} = $j; 
	}
	if ($outType eq 'outTbl') {
		print STDOUT join("\t", @$hd)."\n"; 
	}

	my %h; 
	for my $tl (@{$lines_aref}) {
		my @ta = @$tl; 
		$ta[ $coln{'seq_from'} ] > $ta[ $coln{'seq_to'} ] and @ta[ $coln{'seq_from'}, $coln{'seq_to'} ] = @ta[ $coln{'seq_to'}, $coln{'seq_from'} ]; 
		$opts{'noTrunc'} and do { $ta[ $coln{'trunc'} ] =~ m/^no$/i or next; }; 
		if ($outType eq 'outTbl') {
			print STDOUT join("\t", @ta)."\n"; 
		} elsif ($outType eq 'outBed') {
			print STDOUT join("\t", $ta[ $coln{'target_name'} ], $ta[ $coln{'seq_from'} ]-1, $ta[ $coln{'seq_to'} ], '.', 0, $ta[ $coln{'strand'} ])."\n"; 
		} elsif ($outType eq 'outGff') {
			my $s1 = 'p'; $ta[ $coln{'strand'} ] eq '-' and $s1 = 'm'; 
			my $id = "$ta[ $coln{'target_name'} ]__$ta[ $coln{'seq_from'} ]__$ta[ $coln{'seq_to'} ]__$s1"; 
			my $new_id = $id; 
			my $cnt = 0; 
			while (defined $h{$new_id}) {
				$cnt++; 
				$new_id = "${id}_n$cnt"; 
				$cnt >= 1e6 and &stopErr("[Err] Something maybe wrong for $new_id. \n"); 
			}
			$id = $new_id; 
			print STDOUT join("\t", $ta[ $coln{'target_name'} ], 'infernal', 'rRNA', $ta[$coln{'seq_from'}], $ta[$coln{'seq_to'}], $ta[$coln{'score'}], $ta[$coln{'strand'}], '.', "ID=$id; bias=$ta[$coln{'bias'}]; eval=$ta[$coln{'E-value'}]; inc=\"$ta[$coln{'inc'}]\"; pass=$ta[$coln{'pass'}]; trunc=\"$ta[$coln{'trunc'}]\"; query=\"$ta[$coln{'query_name'}]\"")."\n"; 
		}
	}
}

# --fmt 2 --tblout .tbl
# #idx target name   accession query name           accession clan name mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc olp an
# #--- ------------- --------- -------------------- --------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- --- --
# 1    MIR394        RF00688   Cla97.STRG.659.2     -         -          cm        1      111      313      469      +    no    1 0.32  26.2   61.1   1.1e-11  !   *
# 1    tRNA          RF00005   Cla97.STRG.36651.1   -         CL00001    cm        1       71     1960     1889      -    no    1 0.56   0.0   58.1   9.9e-12  !   *
# 2    tRNA          RF00005   Cla97.STRG.36651.1   -         CL00001    cm        1       71     4231     4160      -    no    1 0.56   0.0   58.1   9.9e-12  !   *
# 3    tRNA          RF00005   Cla97.STRG.36651.1   -         CL00001    cm        1       71      730      659      -    no    1 0.57   0.0   56.7   2.5e-11  !   *
# 4    tRNA          RF00005   Cla97.STRG.36651.1   -         CL00001    cm        1       71     3605     3534      -    no    1 0.54   0.0   55.4   5.7e-11  !   *
# 5    tRNA          RF00005   Cla97.STRG.36651.1   -         CL00001    cm        1       71     1393     1322      -    no    1 0.53   0.0   54.4   1.1e-10  !   *
# 6    tRNA          RF00005   Cla97.STRG.36651.1   -         CL00001    cm        1       71     3148     3077      -    no    1 0.47   0.0   42.4   2.5e-07  !   *
# #
# # Program:         cmscan

#target_name  accession  query_name  accession  mdl  mdl_from  mdl_to  seq_from  seq_to  strand trunc pass gc bias score E-value inc description_of_target
##target name  accession query name           accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of
##------------ --------- -------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- ---------------
#Cmo_Scf00021  -         5S_rRNA              RF00001    cm        1      119  4452908  4453026      +    no    1 0.55   0.0  126.3     7e-27 !   -
sub load_cmscan_tbl_1 {
	my $ifh; 
	if (!(-t)) {
		$ifh = \*STDIN; 
	} elsif (@ARGV >= 1) {
		$ifh = &openFH($ARGV[0], '<'); 
	}
	my @lines; 
	while (<$ifh>) {
		chomp; 
		push(@lines, $_); 
	}
	close($ifh); 
	$lines[1] =~ m!^#\-! or die "The input file format is wrong. line_1: $lines[1]\n"; 
	my @ele_len; 
	{
		my $t1 = $lines[1]; 
		while ($t1 =~ s!^(#?\-+\s*)!!) {
			my $l1 = length($1); 
			push(@ele_len, $l1); 
		}
		$t1 eq '' or die "bad line_1: $lines[1]\n"; 
	}
	my @back_lines; 
	for (my $i=0; $i<@lines; $i++) {
		$i == 1 and next; 
		my $tline = $lines[$i]; 
		if ( $tline =~ m!^#idx|^#target! ) {
			my @ta; 
			for (my $j=0; $j<@ele_len; $j++) {
				my $e1 = substr($tline, 0, $ele_len[$j]); 
				substr($tline, 0, $ele_len[$j]) = ''; 
				$e1 =~ s!^\s*|\s*$!!g; 
				push(@ta, $e1); 
			}
			$ta[0] =~ s!^#+!!; 
			$ta[-1] .= $tline; 
			for my $tb (@ta) {
				$tb =~ s!\s+!_!g; 
			}
			push(@back_lines, [@ta]); 
			next; 
		} elsif ( $tline =~ m!^\s*(#|$)! ) {
			next; 
		} else {
			my @ta; 
			for (my $j=0; $j<@ele_len; $j++) {
				my $e1 = substr($tline, 0, $ele_len[$j]); 
				substr($tline, 0, $ele_len[$j]) = ''; 
				$e1 =~ s!^\s*|\s*$!!g; 
				push(@ta, $e1); 
			}
			$ta[-1] .= $tline; 
			push(@back_lines, [@ta]); 
		}
	}
	return(@back_lines); 
}# load_cmscan_tbl_1() 


