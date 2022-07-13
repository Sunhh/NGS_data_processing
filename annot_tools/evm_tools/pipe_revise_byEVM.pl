#!/usr/bin/perl
# [01/15/2022] Currently I require the gene IDs are unique across all input genomes!!!
# [01/15/2022] Use the information provided by liftoff in gff3 files for filtering.
# [01/19/2022] Adding "-polish" parameter to liftoff will cause this program fail when running, so I prefer to revise the CDS regions later by myself.
#   It is clear that although we have used "-exclude_partial" parameter and applied "valid_ORFs" tags to filter out incomplete genes, there are still gene models with incorrect CDS, in which the stop codons were found in front of tailing bases, or the CDS length can not be divided exactly by 3. 
# [01/20/2022] Instead of using "-polish", add a perl script to fix the stop codon problem.
# [5/25/2022] Rewrite this script to make it simple to read. In this script, I don't implement any function, and I only care about the input parameters.
# [5/26/2022] Select gene models with EVM software, and then add the missing genes by "longCDS" strategy.
# [5/27/2022] It is annoying to have to wait for other liftoff runs finishing, so I want to copy the gff and genome files into my own directory to be independent.

use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;
use ConfigSunhh;
my $cs_obj = ConfigSunhh->new();
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "in_cfg:s", # configuration file
  "in_refList:s",  # 
  "in_addList:s",  # 
  "in_parList:s",  # A file with columns as: tag1 \t tag2 \t identity (-s) \t coverage (-a) \t (-sc)
  "out_pref:s",    # 
  "max_olapLen:i", # 0
  "max_olapRat:f", # 0
  "bestModel:s",   # longCDS
  "keepTmp!",
  "help!",
);

my $htxt = <<HHH;
####################################################################################################
# perl $0 -in_cfg param.cfg -in_refList ref_list  -in_addList add_list  -in_parList param_list  -out_pref merged
# 
# Version 4.2.
# 
#   -bestModel         [longCDS] or EVM
#
#   -max_olapLen       [in_cfg or 0]
#   -max_olapRat       [in_cfg or 0]
# 
#   -keepTmp           [Boolean]
####################################################################################################
HHH

for my $k1 (qw/in_cfg in_refList in_addList in_parList out_pref/) {
  defined $opts{$k1} or &LogInforSunhh::usage($htxt);
}

my %par;
$cs_obj->getConfig( 'cfg_file' => $opts{'in_cfg'}, 'replace' => 1, 'hash_r'=> \%par );
for my $k1 (keys %opts) { $par{$k1} = $opts{$k1}; }
$par{'max_olapLen'} //= 0;
$par{'max_olapRat'} //= 0;
$par{'bestModel'}   //= 'longCDS';

# !!! I should set default values for all parameters in the param.cfg file later.

my $wdir = &fileSunhh::new_tmp_dir('create'=>1);

# Section: Load ref/add lists.
my ($rTag, $fn_rGenom, $fn_rGff, $fn_rPep, $rCN); # Pep is not used yet, but I record it here.
### In liftoff, "r" ("ref") is the target genome onto which the genes will be mapped, while "q" (qry) is source genome providing gene models.
### $rTag        : the tag of the target genome.
### $fn_rGenom   : the fasta file of the target genome.
### $fn_rGff     : the gff-format original gene models of the target genome.
### $fn_rPep     : the fasta file of the target genome original protein sequences. Not used yet.
### $rCN         : the column number of the target genome in a CHR mapping file. Dropped in this script.
my @fn_add;   # This stores information of source genomes.
my %par_list; # {tag1}{tag2} = [ identity (-s), coverage (-a), (-sc) ];
{
  # Load the ref_list. Only the first line is considered.
  my @fn_ref = &fileSunhh::load_tabFile($par{'in_refList'});
  ($rTag, $fn_rGenom, $fn_rGff, $fn_rPep, $rCN) = @{$fn_ref[0]};
  ### Copy input files to my directory.
  &fileSunhh::_copy($fn_rGenom, "$wdir/ref.genome.fa"); $fn_rGenom = "$wdir/ref.genome.fa";
  &fileSunhh::_copy($fn_rGff,   "$wdir/ref.ann.gff3");  $fn_rGff   = "$wdir/ref.ann.gff3";
  -e $fn_rPep and do { &fileSunhh::_copy($fn_rPep,   "$wdir/ref.p.fa");      $fn_rPep   = "$wdir/ref.p.fa"; };
  # Load the add_list (for qry). Each line stands for a source genome and its gene model set.
  @fn_add = &fileSunhh::load_tabFile($par{'in_addList'});
  for (my $i=0; $i<@fn_add; $i++) {
    &fileSunhh::_copy($fn_add[$i][1], "$wdir/add.$i.genome.fa"); $fn_add[$i][1] = "$wdir/add.$i.genome.fa";
    &fileSunhh::_copy($fn_add[$i][2], "$wdir/add.$i.ann.gff3");  $fn_add[$i][2] = "$wdir/add.$i.ann.gff3";
    -e $fn_add[$i][3] and do { &fileSunhh::_copy($fn_add[$i][3], "$wdir/add.$i.p.fa"); $fn_add[$i][3] = "$wdir/add.$i.p.fa"; };
  }
  my @fn_par = &fileSunhh::load_tabFile($par{'in_parList'});
  for (@fn_par) {
    $_->[4] //= 1.0;
    $par_list{$_->[0]}{$_->[1]} = [$_->[2], $_->[3], $_->[4]]; # {tag1}{tag2} = [ -s, -a, -sc ];
    $par_list{$_->[1]}{$_->[0]} = [$_->[2], $_->[3], $_->[4]];
  }
}

# Section: Run liftoff to get transferred gene models from source genomes (addList).
&fileSunhh::write2file("$wdir/gff_list", "$fn_rGff\n", '>'); # gff_list lists all gff-format gene model sets for overlapping check and model selection.
for (my $i=0; $i<@fn_add; $i++) {
  my ($qTag, $fn_qGenom, $fn_qGff, $fn_qPep, $qCN) = @{$fn_add[$i]};
  defined $par_list{$qTag}{$rTag} or &stopErr("[Err] parm {$qTag}{$rTag}\n");
  my ($par_ident, $par_cov, $par_idCopy) = @{$par_list{$qTag}{$rTag}};
  my $cmd_1 = "$par{'exe_liftoff'} $par{'fixPar_liftoff'} -sc $par_idCopy -s $par_ident -a $par_cov -dir $wdir/intermediate.$i -o $wdir/on.$i.ori.gff3 -u $wdir/unmapped.$i -g $fn_qGff $fn_rGenom $fn_qGenom";
  &runCmd($cmd_1);
  $cmd_1 = "$par{'pl_ret_good_model'} $wdir/on.$i.ori.gff3 > $wdir/on.$i.good.gff3";
  &runCmd($cmd_1);
  &fileSunhh::write2file("$wdir/gff_list", "$wdir/on.$i.good.gff3\n", '>>');
}

# Section: Find the best gene model if there are overlapping ones.
#   There are two choices: (1) Select the longest one; (2) Select by EVM with additional information.
### Option (1) Select the longest one.
if ($par{'bestModel'} =~ m!^\s*longCDS\s*$!i) {
  my $cmd_1 = "$par{'pl_get_longCDS_in_gff'} $par{'max_olapLen'} $par{'max_olapRat'} $wdir/gff_list | $par{'pl_deal_table'} -column 3 > $wdir/slct_best.mrnaID";
  &runCmd($cmd_1); $cmd_1 = "";
  my $gff3_keep = join(" ", map { $_->[0] } &fileSunhh::load_tabFile("$wdir/gff_list") );
  &runCmd("cat $gff3_keep | $par{'pl_deal_gff3'} -gffret $wdir/slct_best.mrnaID -idType 'mRNA' > $wdir/slct_best.gff3");
  &runCmd("$par{'pl_rm_HPartFrame'} $wdir/slct_best.gff3 > $par{'out_pref'}.gff3");
  &runCmd("$par{'pl_cds_from_gff3'} -genome_fas $fn_rGenom  -genome_gff $par{'out_pref'}.gff3 > $par{'out_pref'}.c.fa");
  &runCmd("$par{'pl_deal_fasta'} $par{'out_pref'}.c.fa -cds2aa -infer_frame > $par{'out_pref'}.p.fa");
} elsif ($par{'bestModel'} =~ m!^\s*evm\s*$!i) {
  my $oriDir = &fileSunhh::_abs_path("./");
  mkdir("$wdir/evm/");
  my @fn_keepGff = map { $_->[0] } &fileSunhh::load_tabFile("$wdir/gff_list");

  # Create gene_predictions input.
  ### The first file is the current prediction.
  open O1,'>',"$wdir/evm/gene_predictions.comb.gff3" or die "$!\n";
  open G1,'<',"$fn_keepGff[0]" or die "$!\n";
  while (<G1>) { m!^\s*(#|$)! and next; chomp; my @ta=split(/\t/, $_); $ta[1] = "origin"; print O1 join("\t", @ta)."\n"; }
  close G1;
  ### The following gff files are the lift-over predictions.
  for my $fn_t1 (@fn_keepGff[1 .. $#fn_keepGff]) {
    open G2,'<',"$fn_t1" or die "$!\n";
    while (<G2>) { m!^\s*(#|$)! and next; chomp; my @ta=split(/\t/, $_); $ta[1] = "Liftoff"; print O1 join("\t", @ta)."\n"; }
    close G2;
  }
  close O1;
  &runCmd("$par{'pl_fitEvmPred'} $wdir/evm/gene_predictions.comb.gff3 > $wdir/evm/gene_predictions.evm.gff3");

  # Create EVM weight matrix file.
  &fileSunhh::write2file( "$wdir/evm/evm_weight.txt", join("\t", qw/ABINITIO_PREDICTION origin/, $par{'EVM_abPred_weight'})."\n", '>' );
  &fileSunhh::write2file( "$wdir/evm/evm_weight.txt", join("\t", qw/ABINITIO_PREDICTION Liftoff/, $par{'EVM_abPred_weight'})."\n", '>>' );
  for my $l1 (&fileSunhh::load_tabFile($par{'fn_evm_weight'})) {
    $l1->[0] =~ m/PREDICTION/ and next;
    &fileSunhh::write2file( "$wdir/evm/evm_weight.txt", join("\t", @$l1)."\n", '>>' );
  }

  # Copy protein alignment and transcript alignment files and genome file.
  &fileSunhh::_copy($par{'fn_evm_protAln'}, "$wdir/evm/prot_aln.evm.gff3");
  &fileSunhh::_copy($par{'fn_evm_estAln'}, "$wdir/evm/est_aln.evm.gff3");
  &fileSunhh::_copy($fn_rGenom, "$wdir/evm/genome.fa");

  # Run EVM;
  chdir("$wdir/evm/");
  &runCmd("$par{'EVM_HOME'}/EvmUtils/partition_EVM_inputs.pl  --genome genome.fa --gene_predictions gene_predictions.evm.gff3  --protein_alignments prot_aln.evm.gff3  --transcript_alignments est_aln.evm.gff3  $par{'para_evm_partition'} --partition_listing partitions_list.out");
  my $cmd_2 = "$par{'EVM_HOME'}/EvmUtils/write_EVM_commands.pl  --genome genome.fa  --weights `pwd`/evm_weight.txt ";
  $cmd_2 .= " --gene_predictions  gene_predictions.evm.gff3  --protein_alignments prot_aln.evm.gff3  --transcript_alignments est_aln.evm.gff3 ";
  $cmd_2 .= " --output_file_name evm.out  ";
  $cmd_2 .= " --partitions partitions_list.out > commands.list ";
  &runCmd($cmd_2);
  &runCmd("$par{'pl_run_cmd_inBatch'} commands.list -cpuN $par{'cpuN'}");
  &runCmd("$par{'EVM_HOME'}/EvmUtils/recombine_EVM_partial_outputs.pl  --partitions partitions_list.out  --output_file_name evm.out");
  &runCmd("$par{'EVM_HOME'}/EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out  --output evm.out  --genome genome.fa");
  &runCmd("find . -regex \".*evm.out.gff3\" -exec cat \{\} \\; > EVM.all.gff3");
  &runCmd("$par{'pl_rm_HPartFrame'} EVM.all.gff3 > EVM.fixFrame.gff3");
  chdir($oriDir);

  # It seems that EVM will drop gene models for some reason I don't know, so I'd like to draw these dropped genes back.
  mkdir("$wdir/after_evm/");
  &fileSunhh::write2file("$wdir/after_evm/gff_list", '', '>');
  for (my $i=0; $i<@fn_keepGff; $i++) {
    &runCmd("$par{'pl_deal_gff3'} -compare2gffC $wdir/evm/EVM.fixFrame.gff3  -inGff $fn_keepGff[$i] -rmOvlap -rmOvlapLen 1 -rmOvlapType CDS -rmOvlapStrand Single > $wdir/after_evm/novel2evm.$i.gff3");
    &fileSunhh::write2file("$wdir/after_evm/gff_list", "$wdir/after_evm/novel2evm.$i.gff3\n", '>>');
  }
  &runCmd("$par{'pl_get_longCDS_in_gff'} $par{'max_olapLen'} $par{'max_olapRat'} $wdir/after_evm/gff_list | $par{'pl_deal_table'} -column 3 > $wdir/after_evm/slct_best.mrnaID");
  my $gff3_keep = join(" ", map { $_->[0] } &fileSunhh::load_tabFile("$wdir/after_evm/gff_list") );
  &runCmd("cat $gff3_keep | $par{'pl_deal_gff3'} -gffret $wdir/after_evm/slct_best.mrnaID -idType 'mRNA' > $wdir/after_evm/slct_best.gff3");
  &runCmd("$par{'pl_rm_HPartFrame'} $wdir/after_evm/slct_best.gff3 > $wdir/after_evm/slct_best.fixFrame.gff3");

  # Collect results.
  &runCmd("cat $wdir/evm/EVM.fixFrame.gff3  $wdir/after_evm/slct_best.fixFrame.gff3 > $par{'out_pref'}.gff3");
  &runCmd("$par{'pl_cds_from_gff3'} -genome_fas $fn_rGenom -genome_gff $par{'out_pref'}.gff3 > $par{'out_pref'}.c.fa");
  &runCmd("$par{'pl_deal_fasta'} -cds2aa -infer_frame $par{'out_pref'}.c.fa > $par{'out_pref'}.p.fa");
} else {
  &stopErr("[Err] -bestModel [$par{'bestModel'}] not supported yet!\n");
}

# Section: Clean data.
$par{'keepTmp'} or &fileSunhh::_rmtree($wdir);

