#!/usr/bin/perl
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;

!@ARGV and die "perl $0 in_transcript.fa pasa_dbName compreh_prefix\n";

my $inTran = shift;
my $pasaDb = shift;
my $compre = shift;

# Part 1: Fixed file names.
my @p1_rm = qw/blat.spliced_alignments.gff3 gmap.spliced_alignments.gff3 pasa_run.log.dir gmap.spliced_alignments.gff3 blat.spliced_alignments.gff3 gmap.spliced_alignments.gff3.completed blat.spliced_alignments.gff3/;
push(@p1_rm, qw/pblat_outdir cleaning_1 outparts_cln.sort/);
push(@p1_rm, "11.ooc");
for my $a1 (@p1_rm) {
  -e $a1 and &fileSunhh::_rmtree($a1);
}
-e 'alignment.validations.output' and &runCmd("gzip alignment.validations.output");
&runCmd("rm -f blat.spliced_alignments.gff3");

# Part 2: in_transcript.fa related
my @p2_suff = qw/clean.transdecoder_dir clean.transdecoder_dir.__checkpoints clean.transdecoder_dir.__checkpoints_longorfs/;
push(@p2_suff, qw/cidx clean.cidx clean.fai cln clean/);
push(@p2_suff, qw/clean.transdecoder.gff3.fl_accs clean.transdecoder.cds clean.transdecoder.pep clean.transdecoder.bed clean.transdecoder.gff3/);
for my $a2 (@p2_suff) {
  -e "$inTran.$a2" and &fileSunhh::_rmtree("$inTran.$a2");
}
-e $inTran and &fileSunhh::_rmtree($inTran);
-e "err_seqcl_${inTran}.log" and &fileSunhh::_rmtree("err_seqcl_${inTran}.log");
-e "seqcl_${inTran}.log" and &fileSunhh::_rmtree("seqcl_${inTran}.log");

# Part 3: pasa assembly related.
my $pasaAsmFa = "$pasaDb.assemblies.fasta";
-e $pasaAsmFa and &runCmd("gzip $pasaAsmFa");
if (-e "$pasaAsmFa.transdecoder.genome.ID.gff3") {
  &runCmd("gzip $pasaAsmFa.transdecoder.genome.ID.gff3");
  -e "$pasaAsmFa.transdecoder.genome.gff3" and &fileSunhh::_rmtree("$pasaAsmFa.transdecoder.genome.gff3");
} else {
  -e "$pasaAsmFa.transdecoder.genome.gff3" and &runCmd("gzip $pasaAsmFa.transdecoder.genome.gff3");
}
for my $a3 (qw/transdecoder.gff3 transdecoder.genome.ID.gff3 transdecoder.genome.gff3/) {
  -e "$pasaAsmFa.$a3" and &runCmd("gzip $pasaAsmFa.$a3");
}
for my $a3 (qw/transdecoder_dir transdecoder_dir.__checkpoints transdecoder_dir.__checkpoints_longorfs transdecoder.bed transdecoder.cds transdecoder.pep transdecoder.genome.bed transdecoder.genome.ID.gff3.completeID transdecoder.genome.ID.gff3.completeID.p.fa/) {
  -e "$pasaAsmFa.$a3" and &fileSunhh::_rmtree("$pasaAsmFa.$a3");
}
-e "$pasaDb.pasa_assemblies_described.txt" and &runCmd("gzip $pasaDb.pasa_assemblies_described.txt");
-e "$pasaDb.polyAsites.fasta" and &fileSunhh::_rmtree("$pasaDb.polyAsites.fasta");
-e "$pasaDb.pasa_alignment_assembly_building.ascii_illustrations.out" and &fileSunhh::_rmtree("$pasaDb.pasa_alignment_assembly_building.ascii_illustrations.out");
-e "__pasa_${pasaDb}_mysql_chkpts" and &fileSunhh::_rmtree("__pasa_${pasaDb}_mysql_chkpts");
for my $s3 (qw/gtf gff3 bed/) {
  -e "$pasaDb.pasa_assemblies.$s3" and &fileSunhh::_rmtree("$pasaDb.pasa_assemblies.$s3");
  for my $soft (qw/blat gmap/) {
    -e "$pasaDb.valid_${soft}_alignments.$s3" and &fileSunhh::_rmtree("$pasaDb.valid_${soft}_alignments.$s3");
    -e "$pasaDb.failed_${soft}_alignments.$s3" and &fileSunhh::_rmtree("$pasaDb.failed_${soft}_alignments.$s3");
  }
}

# Part 4: comprehensive dataset.
if (-e "$compre/${compre}.4maker.gff3") {
  &runCmd("gzip $compre/${compre}.4maker.gff3");
  &fileSunhh::_rmtree("$compre/${compre}.gff3");
} else {
  -e "$compre/${compre}.gff3" and &runCmd("gzip $compre/${compre}.gff3");
}
for my $a4 (qw/details geneToTrans_mapping/) {
  -e "$compre/${compre}.$a4" and &runCmd("gzip $compre/${compre}.$a4");
}
for my $a4 (qw/bed fasta/) {
  -e "$compre/${compre}.$a4" and &fileSunhh::_rmtree("$compre/${compre}.$a4");
}

