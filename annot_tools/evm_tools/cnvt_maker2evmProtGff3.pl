#!/usr/bin/perl
# [5/24/2022] Convert gff3 files from maker input format to EVM input format.
#   There will be no Parent tag and the ID tag should be the same as the parent ID for all CDS/exon features from the same parent.
#   An example: 
#     Contig1 nap-nr_minus_rice.fasta nucleotide_to_protein_match     8392    8470    50.00   -       .       ID=match.nap.nr_minus_rice.fasta.37;Target=RF|YP_440341.1|83716234|NC_007650 196 222
#     Contig1 nap-nr_minus_rice.fasta nucleotide_to_protein_match     7650    7786    26.09   -       .       ID=match.nap.nr_minus_rice.fasta.37;Target=RF|YP_440341.1|83716234|NC_007650 222 268
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "in_featType:s", # CDS
  "out_featType:s", # nucleotide_to_protein_match
  "out_featSource:s", # spliced_protein_alignments
  "help!",
);
my $htxt = <<HHH;
####################################################################################################
# perl $0  ./src_data/protein/*.spaln.s2.4maker.gff3 > prot_aln.evm.gff3
#
# -in_featType     [CDS,match_part] or exon.
# -out_featType    [nucleotide_to_protein_match]
# -out_featSource  [spliced_protein_alignments] '' means not to change it.
#
HHH

-t and !@ARGV and die "$htxt\n";

$opts{'in_featType'} //= 'CDS,match_part';
$opts{'in_featType'} = lc($opts{'in_featType'});
my %goodFeat = map { s!\s!!g; $_=>1; } split(/,/, $opts{'in_featType'});
$opts{'out_featType'} //= 'nucleotide_to_protein_match';
$opts{'out_featSource'} //= 'spliced_protein_alignments';

while (<>) {
  m!^\s*(#|$)! and next;
  m!^\s*\>! and last;
  chomp;
  my @ta=split(/\t/, $_);
  defined $goodFeat{lc($ta[2])} or next;
  $ta[8] =~ s!^(?i:ID=[^;\s]+;)?\s*Parent=([^;\s]+);!ID=$1;! or die "[Err] Bad line 1: $_\n";
  $ta[2] = $opts{'out_featType'};
  $opts{'out_featSource'} ne '' and $ta[1] = $opts{'out_featSource'};
  print STDOUT join("\t", @ta)."\n";
}


# C31_Chr18       blastx  protein_match   13277446        13284989        2259    +       .       ID=2.1:mRNA00001;Name=AT1G73950.1;
# C31_Chr18       blastx  match_part      13277446        13277561        216     +       0       ID=2.1:cds00001;Parent=2.1:mRNA00001;Target=AT1G73950.1 1 39 +;
# C31_Chr18       blastx  match_part      13278336        13278426        151     +       1       ID=2.1:cds00002;Parent=2.1:mRNA00001;Target=AT1G73950.1 40 69 +;
# C31_Chr18       blastx  match_part      13278546        13278630        157     +       0       ID=2.1:cds00003;Parent=2.1:mRNA00001;Target=AT1G73950.1 70 97 +;
# C31_Chr18       blastx  match_part      13278846        13278930        130     +       2       ID=2.1:cds00004;Parent=2.1:mRNA00001;Target=AT1G73950.1 98 126 +;
# C31_Chr18       blastx  match_part      13279971        13280040        150     +       1       ID=2.1:cds00005;Parent=2.1:mRNA00001;Target=AT1G73950.1 127 149 +;
# C31_Chr18       blastx  match_part      13280856        13280925        140     +       0       ID=2.1:cds00006;Parent=2.1:mRNA00001;Target=AT1G73950.1 150 172 +;
# C31_Chr18       blastx  match_part      13282649        13282694        117     +       2       ID=2.1:cds00007;Parent=2.1:mRNA00001;Target=AT1G73950.1 173 188 +;
# C31_Chr18       blastx  match_part      13282847        13283096        371     +       1       ID=2.1:cds00008;Parent=2.1:mRNA00001;Target=AT1G73950.1 189 270 +;
# C31_Chr18       blastx  match_part      13283337        13283542        270     +       0       ID=2.1:cds00009;Parent=2.1:mRNA00001;Target=AT1G73950.1 271 339 +;
# C31_Chr18       blastx  match_part      13283700        13283774        138     +       1       ID=2.1:cds00010;Parent=2.1:mRNA00001;Target=AT1G73950.1 340 364 +;
# C31_Chr18       blastx  match_part      13284290        13284359        157     +       1       ID=2.1:cds00011;Parent=2.1:mRNA00001;Target=AT1G73950.1 365 387 +;
# C31_Chr18       blastx  match_part      13284458        13284538        168     +       0       ID=2.1:cds00012;Parent=2.1:mRNA00001;Target=AT1G73950.1 388 414 +;
# C31_Chr18       blastx  match_part      13284725        13284804        175     +       0       ID=2.1:cds00013;Parent=2.1:mRNA00001;Target=AT1G73950.1 415 441 +;
# C31_Chr18       blastx  match_part      13284914        13284989        160     +       1       ID=2.1:cds00014;Parent=2.1:mRNA00001;Target=AT1G73950.1 442 466 +;
# C31_Chr02       blastx  protein_match   13317687        13317953        406     +       .       ID=2.1:mRNA00002;Name=AT1G74660.2;
# C31_Chr02       blastx  match_part      13317687        13317953        389     +       0       ID=2.1:cds00015;Parent=2.1:mRNA00002;Target=AT1G74660.2 1 93 +;
# C31_Chr09       blastx  protein_match   15332219        15333909        1929    +       .       ID=2.1:mRNA00003;Name=AT1G67750.1;
# C31_Chr09       blastx  match_part      15332219        15332343        143     +       0       ID=2.1:cds00018;Parent=2.1:mRNA00003;Target=AT1G67750.1 1 42 +;
# C31_Chr09       blastx  match_part      15332437        15333322        1468    +       1       ID=2.1:cds00019;Parent=2.1:mRNA00003;Target=AT1G67750.1 43 335 +;
# C31_Chr09       blastx  match_part      15333691        15333909        349     +       0       ID=2.1:cds00020;Parent=2.1:mRNA00003;Target=AT1G67750.1 336 408 +;
# C31_Chr02       blastx  protein_match   474122  476515  3629    +       .       ID=2.1:mRNA00004;Name=AT1G05150.1;
# C31_Chr02       blastx  match_part      474122  476515  3629    +       0       ID=2.1:cds00021;Parent=2.1:mRNA00004;Target=AT1G05150.1 1 808 +;
# C31_Chr01       blastx  protein_match   5176239 5176492 187     +       .       ID=2.1:mRNA00005;Name=AT1G70895.1;
# C31_Chr01       blastx  match_part      5176239 5176492 109     +       2       ID=2.1:cds00028;Parent=2.1:mRNA00005;Target=AT1G70895.1 20 99 +;
# C31_Chr01       blastx  protein_match   9279399 9280321 930     -       .       ID=2.1:mRNA00006;Name=AT1G24440.2;

