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
for my $a1 (qw/pblat_outdir pasa_run.log.dir gmap.spliced_alignments.gff3 blat.spliced_alignments.gff3 gmap.spliced_alignments.gff3.completed blat.spliced_alignments.gff3/) {
  -e $a1 and &fileSunhh::_rmtree($a1);
}
-e 'alignment.validations.output' and &runCmd("gzip alignment.validations.output");

# Part 2: in_transcript.fa related
my @p2_suff = qw/clean.transdecoder_dir clean.transdecoder_dir.__checkpoints clean.transdecoder_dir.__checkpoints_longorfs/;
push(@p2_suff, qw/cidx clean.cidx clean.fai cln clean/);
push(@p2_suff, qw/clean.transdecoder.gff3.fl_accs clean.transdecoder.cds clean.transdecoder.pep clean.transdecoder.bed clean.transdecoder.gff3/);
for my $a2 (qw/clean.transdecoder_dir clean.transdecoder_dir.__checkpoints clean.transdecoder_dir.__checkpoints_longorfs/) {
  -e "$inTran.$a2" and &fileSunhh::_rmtree("$inTran.$a2");
}
-e $inTran and &fileSunhh::_rmtree($inTran);

# Part 3: pasa assembly related.
my $pasaAsmFa = "$pasaDb.assemblies.fasta";
-e $pasaAsmFa and &runCmd("gzip $pasaAsmFa");
for my $a3 (qw/transdecoder.gff3 transdecoder.genome.ID.gff3 transdecoder.genome.gff3/) {
  -e "$pasaAsmFa.$a3" and &runCmd("gzip $pasaAsmFa.$a3");
}
for my $a3 (qw/transdecoder_dir transdecoder_dir.__checkpoints transdecoder_dir.__checkpoints_longorfs transdecoder.bed transdecoder.cds transdecoder.pep transdecoder.genome.bed transdecoder.genome.ID.gff3.completeID transdecoder.genome.ID.gff3.completeID.p.fa/) {
  -e "$pasaAsmFa.$a3" and &fileSunhh::_rmtree("$pasaAsmFa.$a3");
}
-e "$pasaDb.pasa_assemblies_described.txt" and &runCmd("gzip $pasaDb.pasa_assemblies_described.txt");
-e "$pasaDb.polyAsites.fasta" and &fileSunhh::_rmtree("$pasaDb.polyAsites.fasta");
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

