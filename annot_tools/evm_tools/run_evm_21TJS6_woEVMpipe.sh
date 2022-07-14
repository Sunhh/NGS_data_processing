### Basic functions.
function exe_cmd {
  echo "[$(date)][CMD] $1"
  eval "$1"
  if [[ $? -eq "0" ]]
  then
    echo "[$(date)][CMD_done] $1"
  else
    echo "[$(date)][CMD_err] $1"
  fi
}

function tsmsg {
  echo "[$(date)][Msg] $1"
}

# Parameter settings.
gnPref="21TJS6"
gnFa="/data/Sunhh/wmhifi/analysis/gene_prediction/in_genomes/21TJS6.fa"
mkPref="21TJS6_FalT1"

fnList="/data/Sunhh/wmhifi/analysis/gene_prediction/07_evm/list.gnPref_mkPref_gnFa_gff"

CONDA_activate_liftoff="source /data/Sunhh/install/anaconda3/bin/activate liftoff"
CONDA_deactivate="conda deactivate"

export EVM_HOME=/data/Sunhh/src/annotation/evm/EVidenceModeler-1.1.1/

PL_trim2cds="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/liftoff_tools/fmt_gff_trim2CDS.pl"
PL_retGoodGff="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/liftoff_tools/ret_good_model_from_liftoff_gff3.pl"
PL_fitEVM="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/evm_tools/fit_evm_pred.pl"
PL_cnvtMk2EVM="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/evm_tools/cnvt_maker2evmProtGff3.pl"
PL_rmHeadPartCDS="perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/liftoff_tools/rm_head_partial_frame_inGff.pl"
PL_getCDSgff="perl /home/Sunhh/tools/github/NGS_data_processing/temp/get_cds_from_gff3.pl"

# Prepare genome fasta file
[[ -e db/ ]] || mkdir db/
exe_cmd "cp -p ${gnFa} ./db/${gnPref}.fa"
[[ -e tmp/ ]] || mkdir tmp/

# Prepare predicted gene models.
### From maker (include abinit with IPR domains)
# /data/Sunhh/cpepo/analysis/gene_prediction/run_C31/06_maker/maker_results/3fal2/iprscan/comb_maker_noFa.gff3
GFF_pred="/data/Sunhh/wmhifi/analysis/gene_prediction/06_maker/maker_results/${mkPref}/iprscan/comb_maker_noFa.gff3"

# Run Liftoff to transfer predictions.
exe_cmd "$CONDA_activate_liftoff"
while IFS=$'\t' read -r -a myArr
do
  if [[ ${myArr[1]} == $mkPref ]];
  then
    continue
  fi
  exe_cmd "$PL_trim2cds  ${myArr[3]} > tmp/${myArr[1]}.trim2cds.gff3"
  exe_cmd "liftoff -polish -copies -sc 0.97 -exclude_partial -s 0.97 -a 0.95 -p 40 -dir tmp/im.${myArr[1]} -o tmp/on.${myArr[1]}.ori.gff3 -u tmp/unmapped.${myArr[1]} -g tmp/${myArr[1]}.trim2cds.gff3 ./db/${gnPref}.fa ${myArr[2]}"
  exe_cmd "$PL_retGoodGff  tmp/on.${myArr[1]}.ori.gff3_polished > tmp/on.${myArr[1]}.gff3"
done < "$fnList"
exe_cmd "$CONDA_deactivate"

# Skip evidences from PASA alignments.
### perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/pasa_gff_to_alnGff.pl -notPasa -addTag 'tdcPrj:' ../../db/from_pasa/C31prj_pasa.assemblies.fasta.transdecoder.genome.ID.gff3 | grep -v ^$ > pasaPrj.tdc.gff3
### perl /home/Sunhh/tools/github/NGS_data_processing/annot_tools/pasa_gff_to_alnGff.pl -notPasa -addTag 'tdcSra:' ../../db/from_pasa/C31sra_pasa.assemblies.fasta.transdecoder.genome.ID.gff3 | grep -v ^$ > pasaSra.tdc.gff3

# Combine gene predictions.
cat $GFF_pred > tmp/gene_predictions.evm.gff3
while IFS=$'\t' read -r -a myArr
do
  if [[ ${myArr[1]} == $mkPref ]];
  then
    continue
  fi
  cat tmp/on.${myArr[1]}.gff3 >> tmp/gene_predictions.evm.gff3
done < "$fnList"
$PL_fitEVM tmp/gene_predictions.evm.gff3 > gene_predictions.evm.gff3
rm -f tmp/gene_predictions.evm.gff3

# Prepare protein alignments.
cat /data/Sunhh/wmhifi/analysis/gene_prediction/03_aln_prot/run_${gnPref}/*.spaln.s2.4maker.gff3 | $PL_cnvtMk2EVM > prot_aln.evm.gff3

# Prepare EST alignments.
### Skip this for now.
# cat ../../db/from_pasa/C31_compreh_prj.4maker.gff3 ../../db/from_pasa/C31_compreh_sra.4maker.gff3 | perl cnvt_maker2evmProtGff3.pl -out_featSource "" -out_featType "EST_match" > est_aln.evm.gff3

# Prepare weight matrix file.
### evm_weight.txt
echo -e "ABINITIO_PREDICTION\tmaker\t6" > evm_weight.txt
echo -e "ABINITIO_PREDICTION\tLiftoff\t5" >> evm_weight.txt
#echo -e "OTHER_PREDICTION\ttransdecoder\t3" >> evm_weight.txt
echo -e "PROTEIN\tspliced_protein_alignments\t1" >> evm_weight.txt
#echo -e "TRANSCRIPT\tblat-DB_pasa\t1" >> evm_weight.txt
#echo -e "TRANSCRIPT\tgmap-DB_pasa\t1" >> evm_weight.txt
#echo -e "TRANSCRIPT\tassembler-DB_pasa" >> evm_weight.txt

# Run EVM
[[ -e evm1 ]] || mkdir evm1/
cp -p evm_weight.txt evm1/

exe_cmd "$EVM_HOME/EvmUtils/partition_EVM_inputs.pl  --genome $gnFa  --gene_predictions gene_predictions.evm.gff3  --protein_alignments prot_aln.evm.gff3  --segmentSize 800000  --overlapSize 10000  --partition_listing evm1/partitions_list.out"
# --transcript_alignments est_aln.evm.gff3

exe_cmd "$EVM_HOME/EvmUtils/write_EVM_commands.pl  --genome $gnFa  --weights  `pwd`/evm1/evm_weight.txt  --gene_predictions gene_predictions.evm.gff3  --protein_alignments prot_aln.evm.gff3  --output_file_name evm.out --partitions evm1/partitions_list.out > evm1/commands.list"
# --transcript_alignments est_aln.evm.gff3

run_cmd_in_batch.pl evm1/commands.list -cpuN 40 > evm1/scrn.commands.list

exe_cmd "$EVM_HOME/EvmUtils/recombine_EVM_partial_outputs.pl --partitions evm1/partitions_list.out  --output_file_name evm1/evm.out"
exe_cmd "$EVM_HOME/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions evm1/partitions_list.out  --output evm1/evm.out  --genome $gnFa"

find . -regex "evm1/.*evm.out.gff3" -exec cat {} \; > evm1/EVM.all.gff3
exe_cmd "$PL_rmHeadPartCDS evm1/EVM.all.gff3 > evm1/EVM.fixFrame.gff3"

exe_cmd "$PL_getCDSgff -genome_fas $gnFa  -genome_gff evm1/EVM.fixFrame.gff3  > evm1/EVM.fixFrame.gff3.cds.fa"

exe_cmd "deal_fasta.pl evm1/EVM.fixFrame.gff3.cds.fa -cds2aa -infer_frame > evm1/EVM.fixFrame.gff3.prot.fa"

# Collect results.

